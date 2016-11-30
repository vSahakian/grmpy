####Run the residual analysis####
'''
Run the Residual analysis for a given parameter file
VJS 7/2016
'''


#Initialize the folders
def init(home,run_name):
    '''
    Initialize the folders for output for a new run (new database/model combo)
    Input:
        home:           String to path of the home environment
        run_name:       String with the name for this db/model combo
    Output: 
        Makes directories...
    '''
    
    from shutil import rmtree,copy
    from os import path,makedirs,environ
    
    clob='y'
    run_dir=path.expanduser(home+run_name+'/')
    if path.exists(run_dir):    #If the path exists...
        print 'Run directory is: '+run_dir
        clob=raw_input('The run directory already exists, are you sure you want to clobber it? (y/n)')
        if clob is 'y' or clob is 'Y':
            clob=raw_input('This deletes everything in '+run_dir+', are you sure you want to do that? (y/n)')
            if clob is 'y' or clob is 'Y':
                rmtree(run_dir)
            else:
                #Don't do anything
                print 'Aborting mission...not clobbering...'
        else:
            #AGain, don't do anything
            print 'Aborting mission...not clobbering...'
    if clob is 'y' or clob is 'Y':
        #Make the main directory
        makedirs(run_dir)
        #And the subdirectories:
        makedirs(run_dir+'event_objs/')
        makedirs(run_dir+'sta_objs/')
        makedirs(run_dir+'E_resids/')
        makedirs(run_dir+'W_resids/')
        makedirs(run_dir+'site_resids/')
        makedirs(run_dir+'figs/')
        makedirs(run_dir+'figs/pdfs/')
        makedirs(run_dir+'raytracing/')
        
        #If it's ok to go ahead and clobber, set the run variable to 1:
        runall=1
        
    if clob is 'n' or clob is 'N':
        #Do not run!!!!
        #Set the run variable to 0:
        runall=0
        
        
    return runall
        
        
def get_total_res(home,run_name,dbpath,modelpath,Mc,ffdf_flag,resaxlim):
    '''
    Get the total residuals for a database.
    Input:
        home:            STring with home path
        run_name:        String with name of the run
        dbpath:          String with path to a pickle database object
        modelpath:       String with path to a gmpe model pickle object
        Mc:              Value aroudn which to center the magnitude squred term (like 8.1 or 8.5, etc.)
        ffdf_flag:       Flag for mag dependent ffdf.  0=off, 1=on
        resaxlim:        Array with axis limits for the resid plot: [[xmin,xmax],[ymin,ymax]]
        
    Output:
        total_residuals:    Array with total residuals
        mw:                 Array, same dim as total_residuals, with matching mw
        mean_residual:      Mean residual
        std_dev_resid:      Standard deviation of the residuals
    '''
    
    import pickle
    import gmpe as gm
    import rescomp as rcomp
    from numpy import where,array,mean,log,log10
    from os import path
    import cdefs as cdf
    from matplotlib.pyplot import savefig
    
    run_dir=path.expanduser(home+run_name+'/')
    
    
    ##########Open the database object:###################
    #Filename:
    dname=dbpath
    datobj=open(dname,'r')
    db=pickle.load(datobj)
    datobj.close()
    
    ##Residual Computation####
    #Load in the model to use:
    mname=modelpath
    datobj=open(mname,'r')
    model=pickle.load(datobj)
    datobj.close()
    
    print 'Opened model %s' % modelpath
    
    #Overall residual, 
    #In some places, vs30 is 0.  Set these to vref.
    vref=760
    mdep_ffdf=ffdf_flag
    #Where are they 0?
    vs30_0ind=where(db.vs30==0)[0]
    #Get vs30 from database...
    vs30=db.vs30
    #Set 0 entries to vref:
    vs30[vs30_0ind]=vref
    
    
    #####Get Total Residuals######
    
    #Now compute the predicted value of PGA...
    d_predicted=gm.compute_model(model.m,model.rng,db.mw,db.r,db.ffdf,vs30,Mc,vref,mdep_ffdf)
    
    # Test the mean from the inversion model:
    #testresiduals = model.G.dot(model.m) - model.d
    testresiduals = d_predicted - model.d
    meantest = mean(testresiduals)
    print 'Mean residual from G and d matrices are %f' % meantest
    
    #Test that the correted data, outside of rcomp, is the same as what is from inversion:
    d_comparison = model.d - (log(db.pga_pg) - 0.6*log(vs30/vref))
    print d_comparison
    
    #Get residuals:
    total_residuals,mean_residual,std_dev=rcomp.total_residual(db,d_predicted,vref)
    
    #Plot residuals and save plots:
    allresid=cdf.total_residuals(db.mw,total_residuals,mean_residual,std_dev)
    f1=allresid.plt_resids(run_name,resaxlim)
    
    #Figure names:
    f1name=run_dir+'figs/'+run_name+'_total.png'
    f1pdf=run_dir+'figs/pdfs/'+run_name+'_total.pdf'

    
    #Save:
    f1.savefig(f1name)
    f1.savefig(f1pdf)

    print 'Saved figures to %s and %s' % (f1name,f1pdf)

    return db.mw,allresid,mean_residual,std_dev
    
    
    
    
def getEW_makeEvents(home,run_name,dbpath,modelpath,Mc,vref,ffdf_flag,resaxlim):
    '''
    Get the Event and Within-Event Residuals, and store in Event objects
    Input:
        home:            STring with home path
        run_name:        String with name of the run
        dbpath:          String with path to a pickle database object
        modelpath:       String with path to a gmpe model pickle object
        Mc:              Number around which to center magnitude squred term (i.e., 8.1 or 8.5)
        vref:            Reference vs30 velocity
        ffdf_flag:       Flag for mag dependent ffdf.  0=off, 1=on
        resaxlim:        Array with axis limits for the resid plot: [[xmin,xmax],[ymin,ymax]]
    Output:
        E_evnum:                Array of event numbers
        E_mw:                   Array of event magnitudes
        E_residual:             Array of event residuals for db
        E_mean:                 Mean event residual for db
        E_std_dev:              Standard deviation for all event residuals in db
        Event_list_object:      An object, saved to run_dir/event_objs/run_name.pckl, 
                                with a list of all event objects for the databases
    '''
    
    import numpy as np
    import cdefs as cdf
    import rescomp as rcomp
    from os import path
    import gmpe as gm
    import cPickle as pickle
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator

    
    #Get directory for output:
    run_dir=path.expanduser(home+run_name+'/')
    eo_dir=run_dir+'event_objs/'
    fig_dir=run_dir+'figs/'
     
    ###Get Event and Within-Event Residuals##\
    #Import the database and model predictions
    ##########Open the database object:###################
    #Filename:
    dname=dbpath
    datobj=open(dname,'r')
    db=pickle.load(datobj)
    datobj.close()
    
    ##Residual Computation####
    #Load in the model to use:
    mname=modelpath
    datobj=open(mname,'r')
    model=pickle.load(datobj)
    datobj.close()
    
    #In some places, vs30 is 0.  Set these to vref.
    print 'vref is %i' % vref
    mdep_ffdf=ffdf_flag
    #Where are they 0?
    vs30_0ind=np.where(db.vs30==0)[0]
    #Get vs30 from database...
    vs30=db.vs30
    #Set 0 entries to vref:
    vs30[vs30_0ind]=vref
    
    #Now compute the predicted value of PGA...
    d_predicted_ln=gm.compute_model(model.m,model.rng,db.mw,db.r,db.ffdf,vs30,Mc,vref,mdep_ffdf)

    ###
    ##Get unique events:
    unique_events=np.unique(db.evnum)

    ##Make a class for each event; append them to a list, made empty here:
    event_list=[]
    d_predicted_list=[]
    
    ##Zero out arrays for info used for plotting, for all events...
    E_evnum=[]
    E_mw=[]
    E_residual=[]
    
    ##Zero out lists to use for indexing - use this resulting index, the station index,
    # to sort out the events for the station object later on...
    unique_stnums=np.unique(db.stnum)
    unique_stas=np.unique(db.sta)

    #Loop through the unique events, make each into an object, append to event list
    for i in range(len(unique_events)):
        unique_ind=np.where(unique_events[i]==db.evnum)[0]
        #Get the predicted data for this event only, for each recording
        #d_predicted_i is an array, with len = # of recordings for event):
        d_predicted_i=d_predicted_ln[unique_ind]
    
        #Get the database info for this event:
        evnum_i=db.evnum[unique_ind]
        sta_i=db.sta[unique_ind]
        stnum_i=db.stnum[unique_ind]
        ml_i=db.ml[unique_ind]
        mw_i=db.mw[unique_ind]
        pga_i=db.pga[unique_ind]
        pgv_i=db.pgv[unique_ind]
        pga_pg_i=db.pga_pg[unique_ind]
        r_i=db.r[unique_ind]
        ffdf_i=db.ffdf[unique_ind]
        md_ffdf_i=db.md_ffdf[unique_ind]
        lat_i=db.elat[unique_ind]
        lon_i=db.elon[unique_ind]
        depth_i=db.edepth[unique_ind]
        stlat_i=db.stlat[unique_ind]
        stlon_i=db.stlon[unique_ind]
        stelv_i=db.stelv[unique_ind]
        source_i=db.source_i[unique_ind]
        receiver_i=db.receiver_i[unique_ind]
        

        #The only one that is not directly from teh database is vs30; this is because
        #there were some 0's in it, so above that is compensated for by changing 
        #them to vref.
        vs30_i=vs30[unique_ind]
    
    
        #Make the event object:
        eventi=cdf.event(evnum_i,sta_i,stnum_i,ml_i,mw_i,pga_i,pgv_i,pga_pg_i,r_i,vs30_i,ffdf_i,md_ffdf_i,lat_i,lon_i,depth_i,stlat_i,stlon_i,stelv_i,source_i,receiver_i)
        
        #Get total residual to store:
        total_residual=(np.log(pga_pg_i) - 0.6*np.log(vs30_i/vref))-d_predicted_i
        
        #Add total residual to object:
        eventi.add_total_resid(total_residual)
        
        #Compute the event terms:
        evnum_i,evmw_i,E_residual_i,std_dev_i=rcomp.event_residual(eventi,d_predicted_i)
        
        #Add the residual information to the event object:
        eventi.add_E_resid(E_residual_i,std_dev_i)
        
        #Get the Within-Event Residuals:
        evnum_i,evmw_i,sta_i,stnum_i,W_residuals_i,W_mean_i,W_std_dev_i=rcomp.within_event_residual(eventi,d_predicted_i,eventi.E_residual)    
    
        #Add the within-event residuals to the event object:
        eventi.add_W_resids(W_residuals_i,W_mean_i,W_std_dev_i)
        
        #Append the event object, and the d_predicted to the list:
        event_list.append(eventi)
        d_predicted_list.append(d_predicted_i)
        
        #Append to the residual arrays, to use later for plotting:
        E_evnum.append(evnum_i)
        E_mw.append(evmw_i)
        E_residual.append(E_residual_i)
        
    #Turn those into arrays:
    E_evnum=np.array(E_evnum)
    E_mw=np.array(E_mw)
    E_residual=np.array(E_residual)
    
    #Get stats:
    E_mean=np.mean(E_residual)
    E_std_dev=np.std(E_residual)
    
    #Write the whole event_list list into a pickle object
    fname=eo_dir+run_name+'.pckl'
    flist=open(fname,'w')
    #Loop over each event in event_list, and dump it into the pickle file:
    for i in range(len(event_list)):
        pickle.dump(event_list[i],flist)
    #close the file
    flist.close()
    
    
    ####Plotting####
    
    #Get the file path of the output figure:
    f1figname=fig_dir+run_name+'_E_resids.png'
    f1pdf=fig_dir+'pdfs/'+run_name+'_E_resids.pdf'
    
    f2figname=fig_dir+run_name+'_E_resids_hist.png'
    f2pdf=fig_dir+'pdfs/'+run_name+'_E_resids_hist.pdf'
    
    ##
    #Set up plot to have histogram adjacent to scatter...
    ##
    
    #definitions for axes
    fudge_factor=0.02
    left, width=0.1, 0.70
    bottom, height=0.1,0.81
    left_h=left+width+fudge_factor
    width_h=0.15

    #For titles:
    E_std_dev_short=np.around(E_std_dev,decimals=2)
    
    #Axes:
    axxlim=resaxlim[0]
    axylim=resaxlim[1]

    #Color:
    rgb=np.array([111,168,163])/255.0
    rgb.astype(float)
    color=rgb*np.ones((len(E_mw),3))



    #Plot...    
    f1=plt.figure()
    
    #define axis limits for scatter, and histogram:
    rect_scatter=[left,bottom,width,height]
    rect_histy=[left_h,bottom,width_h,height]
    
    #define axis tick locations for histogram:
    hist_xLocator=MultipleLocator(500) 
    
    #define axes:
    axScatter=plt.axes(rect_scatter)
    axHisty=plt.axes(rect_histy)
    
    axScatter.scatter(E_mw,E_residual,edgecolors=color,facecolors='none',lw=0.8)

    #Histogram:
    #want 4x as many bins as main plot y-axis limit units:
    nbins=(axylim[1]-axylim[0])*4
    
    #set the number of bins, adn the range to be the x axis limits (same as y axis, ln residuals):
    axHisty.hist(E_residual,bins=nbins,range=[axylim[0],axylim[1]],orientation='horizontal',color=rgb)
    
    #Also plot a dashed line at 0:
    axScatter.plot(axxlim,[0,0],linestyle='--',color='0.75')     
    
    #set axis limits:
    #scatter
    axScatter.set_xlim(axxlim)
    axScatter.set_ylim(axylim)
    #histogram
    axHisty.set_ylim(axylim)
    #set axis ticks:
    axHisty.xaxis.set_major_locator(hist_xLocator)
    #set no labels on the y axis:
    axHisty.yaxis.set_ticklabels('')      
    
    #Add titles, limits...
    axScatter.set_xlabel(r"$\mathbf{M}$")
    axScatter.set_ylabel('ln Residual')
    ptitle=r"Event Residuals"+"\n"+"Run: "+run_name+", Std Dev: "+np.str(E_std_dev_short)
    axScatter.set_title(ptitle)
    
    #Show...
    f1.show()
    
    #Save:
    f1.savefig(f1figname)
    f1.savefig(f1pdf)
    
#    #Histogram...
#    f2=plt.figure()
#    #Want 4x as many bins as main plot y-axis limit units:
#    nbins=(axylim[1]-axylim[0])*4
#    #Set the number of bins, adn the range to be the x axis limits (same as y axis, ln residuals):
#    plt.hist(E_residual,bins=nbins,range=[axylim[0],axylim[1]],color=rgb)
#
#    
#    plt.xlabel('ln Residuals')
#    plt.ylabel('Number of occurences')
#            
#    ptitle=r"Event Residuals"+"\n"+"Mean: "+str(np.around(E_mean,decimals=2))+" Std Dev: "+str(np.around(E_std_dev,decimals=2))
#    plt.title(ptitle)
#    
#    f2.show()
#    f2.savefig(f2figname)
#    f2.savefig(f2pdf)
    
    return E_evnum,E_mw,E_residual,E_mean,E_std_dev
    
    
def sta_list(home,run_name,dbfile):
    '''
    Read in a database and an event list object, sort out by station, and 
    computes site residuals, creates station objects, to append to a list.
    Input:
        home:            String with home path
        run_name:        String with name of the run
        dbfile:          String with path to a pickle database object 
    Output: 
        station_list_object:    Writes a list of station objects into an object,
                                into directory /home/run_name/sta_objs/
    '''
    
    import cdefs as cdf
    import numpy as np
    import cPickle as pickle
    import dread
    from os import path
    
    #Get directory for output:
    run_dir=path.expanduser(home+run_name+'/')
    so_dir=run_dir+'sta_objs/'
    
    #Get event list directory for this run:
    eobjfile=run_dir+'event_objs/'+run_name+'.pckl'
    
    #REad database object:
    #Filename:
    fname=dbfile
    datobj=open(fname,'r')
    db=pickle.load(datobj)
    datobj.close()
    
    ##Zero out lists to use for indexing - use this resulting index, the station index,
    # to sort out the events for the station object later on...
    ##Get unique events:
    unique_events=np.unique(db.evnum)
    unique_stnums=np.unique(db.stnum)
    unique_stas=np.unique(db.sta)
    
    #REad event list:
    event_list=dread.read_obj_list(eobjfile)

    ###Make Station Objects, for other residuals####

    #Start an empty list, store all station objects here in the end:    
    station_list=[]
    site_terms=[]
    
    #First loop over the list of unique stations    
    for sta_ind in range(len(unique_stnums)):
        #Station number?
        station_num_i=unique_stnums[sta_ind]
        
        #Station name:
        station_name_i=unique_stas[sta_ind]
        
        #Zero out the lists/arrays that will be used to add to the station object for htis statioN:
        evnum=[]
        ml=[]
        mw=[]
        pga_pg=[]
        pga=[]
        pgv=[]
        r=[]
        ffdf=[]
        md_ffdf=[]
        elat=[]
        elon=[]
        edepth=[]
        stlat=[]
        stlon=[]
        stelv=[]
        source_i=[]
        receiver_i=[]
        #Total residual:
        total_residual=[]
        #Event residuals:
        E_residual=[]
        #Within-event residuals:
        W_residual=[]
        
        #Loop over the unique events:  
        for event_ind in range(len(unique_events)):
            #What event is this?
            eventi=event_list[event_ind]
            #What stations record thsi event?
            event_sta_i=eventi.stnum
            #Is the station being referenced int eh outer loop contained here?
            sta_ev_ind=np.where(event_sta_i==station_num_i)[0]
            
            #If this station records thsi event, grab the information:
            if sta_ev_ind.size!=0:
                vs30=eventi.vs30[sta_ev_ind]
                evnum.append(eventi.evnum[sta_ev_ind])
                ml.append(eventi.ml[sta_ev_ind])
                mw.append(eventi.mw[sta_ev_ind])
                pga_pg.append(eventi.pga_pg[sta_ev_ind])
                pgv.append(eventi.pgv[sta_ev_ind])    
                r.append(eventi.r[sta_ev_ind])    
                ffdf.append(eventi.ffdf[sta_ev_ind])
                md_ffdf.append(eventi.md_ffdf[sta_ev_ind])
                elat.append(eventi.elat[sta_ev_ind])
                elon.append(eventi.elon[sta_ev_ind])
                edepth.append(eventi.edepth[sta_ev_ind])
                stlat.append(eventi.stlat[sta_ev_ind])
                stlon.append(eventi.stlon[sta_ev_ind])
                stelv.append(eventi.stelv[sta_ev_ind])
                source_i.append(eventi.source_i[sta_ev_ind])
                receiver_i.append(eventi.receiver_i[sta_ev_ind])
                total_residual.append(eventi.total_residual[sta_ev_ind])
                E_residual.append(eventi.E_residual)
                W_residual.append(eventi.W_residuals[sta_ev_ind])
                
            elif sta_ev_ind.size==0:
                continue
                
            
        #AFter looping over all events, convert these lists into arrays, to put
        #into a station object:
        #Vs30 stays as just one number
        evnum=np.array(evnum)
        ml=np.array(ml)
        mw=np.array(mw)
        pga_pg=np.array(pga_pg)
        pgv=np.array(pgv)
        r=np.array(r)
        ffdf=np.array(ffdf)
        md_ffdf=np.array(md_ffdf)
        elat=np.array(elat)
        elon=np.array(elon)
        edepth=np.array(edepth)
        stlat=np.array(stlat)
        stlon=np.array(stlon)
        stelv=np.array(stelv)
        source_i=np.array(source_i)
        receiver_i=np.array(receiver_i)
        total_residual=np.array(total_residual)
        E_residual=np.array(E_residual)
        W_residual=np.array(W_residual)
        
        #Put into a station object...
        station_i=cdf.station(station_name_i,station_num_i,vs30[0],evnum,ml,mw,pga_pg,pga,pgv,ffdf,md_ffdf,elat,elon,edepth,stlat,stlon,stelv,source_i,receiver_i,total_residual,E_residual,W_residual)
        
        #Get the site residual:
        station_i.get_site_resid()
        
        #Append this site term and station info to the site_terms list:
        station_i_site_term=station_i.site_resid
        
        
        site_terms.append(station_i_site_term)
        
        #Append to the station list...
        station_list.append(station_i)


    #Write out the station list to an object:
    fname=so_dir+run_name+'.pckl'
    flist=open(fname,'w')
    #Loop over each event in event_list, and dump it into the pickle file:
    for i in range(len(station_list)):
        pickle.dump(station_list[i],flist)
    #close the file
    flist.close()
        
    return station_list



#####
def plot_Wresid(home,run_name,resaxlim):
    '''
    '''
    
    from os import path
    import matplotlib.pyplot as plt
    import numpy as np
    import dread
    
    #Get paths to things, like the object list:
    #Get directory for output:
    run_dir=path.expanduser(home+run_name+'/')
    so_dir=run_dir+'sta_objs/'
    fig_dir=run_dir+'figs/'
    
    #Get the file path of the output figure:
    figname=fig_dir+run_name+'_W_resids.png'
    pdfname=fig_dir+'pdfs/'+run_name+'_W_resids.pdf'
    
    #Get event list directory for this run:
    sobjfile=so_dir+run_name+'.pckl'
    
    #Read in the station list:
    station_list=dread.read_obj_list(sobjfile)
    
    ###Plotting####
    #Plot the event residuals for each station in a different color...
    
    #Axes:
    axxlim=resaxlim[0]
    axylim=resaxlim[1]
    
    #Get the range for the colorbar:
    crangemax=len(station_list)
    crange=np.array(range(len(station_list))).astype(float)/crangemax
    
    #Get the colorbar:
    colors=plt.cm.rainbow(crange)
    
    #Make W_residuals array by appending to use to get std deviation:
    W_res_array=np.array([])
    
    #Plot
    plt.figure()
    for station_ind in range(len(station_list)):
        #Get the station 
        station_i=station_list[station_ind]
        
        #Get what you're plotting...
        mw=station_i.mw
        W_residuals=station_i.W_residual
        
        if station_ind==0:
            W_res_array=W_residuals
        else:
            W_res_array=np.r_[W_res_array,W_residuals]
        
        #Get the color info and label info:
        color_i=colors[station_ind]
        sta_lab=station_i.sta
        
        print color_i
        
        #Plot
        plt.scatter(mw,W_residuals,edgecolors=color_i,facecolors='none',lw=0.8,label=sta_lab)
        plt.hold(True)
     
    #Also plot a line at 0:
    plt.plot(axxlim,[0,0],linestyle='--',color='0.75')
      
    #Get stats:
    W_mean=np.mean(W_res_array)
    
    W_std_dev=np.std(W_res_array)
    W_std_dev_short=np.around(np.std(W_std_dev),decimals=2)
    
    #Add legend:
    plt.legend(loc=4)
    
    #Add titles, limits...
    plt.xlim(axxlim)
    plt.ylim(axylim)
    plt.xlabel(r"$\mathbf{M}$")
    plt.ylabel('ln Residual')
    ptitle=r"Within-Event Residuals by station"+"\n"+"Run: "+run_name+", Std Dev: "+np.str(W_std_dev_short)
    plt.title(ptitle)
    
    #Show...
    plt.show()
    
    #Save:
    plt.savefig(figname)
    plt.savefig(pdfname)
    
    return W_mean,W_std_dev
    
    


#######
def plot_site_WE(home,run_name,resaxlim):
    '''
    Plot the within-event residuals by station
    Input:
        home:               String with path to project home
        run_name:           String with db/model combo
        resaxlim:           Axis limits: [[xmin,xmax],[ymin,ymax]]
    Output:
        
    '''
    
    import cdefs as cdf
    import matplotlib.pyplot as plt
    from os import path
    import dread
    from numpy import ceil,array,floor,remainder
    
    #Get directories:
    run_dir=path.expanduser(home+run_name+'/')
    so_dir=run_dir+'sta_objs/'
    fig_dir=run_dir+'figs/'
    
    #Get the file path of the output figure:
    figname=fig_dir+run_name+'_site_W_resids.png'
    pdfname=fig_dir+'pdfs/'+run_name+'_site_W_resids.pdf'
    
    #Get event list directory for this run:
    sobjfile=so_dir+run_name+'.pckl'
    
    #Read in the station list:
    station_list=dread.read_obj_list(sobjfile)

    #Get info for the plot configuration:
    nsta=len(station_list)
    #Want 3 subplots:
    ncols=3
    #How many rows?
    nrows=ceil(float(nsta)/ncols).astype('int64')
    #Get colors:
    colors=plt.cm.rainbow(array(range(nsta)).astype(float)/nsta)
    
    #Initiate plot:
    fsize=(ncols*4,nrows*2.1)
    fig,axarr=plt.subplots(nrows,ncols,figsize=fsize)    
    
    #Now plot...
    for station_ind in range(len(station_list)):
        #What subplot are we on?
        ax_row=floor(float(station_ind)/ncols).astype('int')
        ax_col=remainder(station_ind,ncols)
        axis_i=axarr[ax_row,ax_col]
        
        #color for plotting:
        color_i=colors[station_ind]
        #axis labeling info:
        if ax_col==0.0:
            ytoggle='on'
        else:
            ytoggle='off'
        
        if ax_row==len(axarr)-1:
            xtoggle='on'
        else:
            xtoggle='off'
        
        #station:
        station_i=station_list[station_ind]
        
        #Now plot:
        station_i.plot_site_WE(axis_i,color_i,resaxlim,xtoggle,ytoggle)
        
        last_row=ax_row
        last_col=ax_col
    #
    ##Disappera the remaining axes unused:
    if last_col<ncols-1:
        #total columnns in last row should be the total number of columns, ncol.
        #How many columns leftover?
        remaining_cols=(ncols-1)-last_col
    
        for col in range(remaining_cols):
            axarr[last_row,((ncols-1)-col)].axis('off')
    
    fig.show()
    
    #Save the figure:
    fig.savefig(figname)
    fig.savefig(pdfname)
    

#######
def get_path_resid_make_object(home,run_name,dbpath,axlims,axlims_dist):
    '''
    Get the path term, meanwhile make an all residuals object
    VJS 8/2016
    Input: 
        home:       Path to the home directory
        run_name:   STring with the run name
        dbpath:     Path to the database being used
        axlims:     Axis limits for plotting [[xmin,xmax],[ymin,ymax]]
        axlims_dist: Axis limits for plotting the distance path terms [[xmin,xmax],[ymin,ymax]]
    Output:
        f_mw:           Figure of path term vs. mw
        f_dist:         Figure of path term vs. r
        allresiduals:   An object containing all db info, all residuals, and raytracing indices
        pterm_mean:     Mean of path term
        pterm_std:      Standard deviation of path term
    '''
    
    from os import path
    import cdefs as cdf
    import cPickle as pickle
    from numpy import mean,std
    
    #Get the run directory:
    run_dir=path.expanduser(home+run_name+'/')
    
    #Neet the directories to the event list, and station list (eo, so):
    eo_dir=run_dir+'event_objs/'
    so_dir=run_dir+'sta_objs/'
    fig_dir=run_dir+'figs/'
    pdf_dir=fig_dir+'pdfs/'

    #Figure names:
    mwfile=fig_dir+run_name+'_pathterms_mw.png'
    mwpdf=pdf_dir+run_name+'_pathterms_mw.pdf'
    
    distfile=fig_dir+run_name+'_pathterms_r.png'
    distpdf=pdf_dir+run_name+'_pathterms_r.pdf'

    #Now get the path for these lists of objects:
    event_list_path=eo_dir+run_name+'.pckl'
    station_list_path=so_dir+run_name+'.pckl'
    
    #Get residuals object:
    allresiduals=cdf.residuals(dbpath,event_list_path,station_list_path)
    
    #Save the residuals object:
    rfilename=run_dir+run_name+'_robj.pckl'
    rfile=open(rfilename,'w')
    pickle.dump(allresiduals,rfile)
    rfile.close()
    
    #Get path term stats:
    pterm_mean=mean(allresiduals.path_terms)
    pterm_std=std(allresiduals.path_terms)
    
    ##Plot:
    #First by magnitude:
    f_mw=allresiduals.plot_path_term_mw(run_name,axlims)
    #save:
    f_mw.savefig(mwfile)
    f_mw.savefig(mwpdf)
    
    #then by distance:
    f_dist=allresiduals.plot_path_term_r(run_name,axlims_dist)
    #save
    f_dist.savefig(distfile)
    f_dist.savefig(distpdf)
    
    #REturn the allr esiduals object:
    return f_mw,f_dist,allresiduals,pterm_mean,pterm_std
 
######
def write_stats(home,run_name,mean_tot,std_dev_tot,E_mean,E_std_dev,W_mean,W_std_dev,pterm_mean,pterm_std):
    '''
    Write out statistics to a file
    Input:
        home:           String to path where the home dir is
        run_name:       String with the name of this run
        mean_tot:       Mean Total residual
        std_dev_tot:    Std Dev Total residual
        E_mean:         Mean Event Residual
        E_std_dev:      Std Dev Event Residual
        W_mean:         Mean Within-Event Residual
        W_std_dev:      Std Dev Within-Event Residual
        pterm_mean:     Mean of the path terms
        pterm_std:      Std Dev of the path terms
    Output:
        STatistics file to home/run_name/run_name_stats.txt 
    '''
    from os import path
    from numpy import str,unique,array,around,mean,std
    import cPickle as pickle
    
    #Make stats file:   
    statsdir=path.expanduser(home+run_name+'/')
    statsfile=statsdir+run_name+'_stats.txt'
    
    #Open the residuals object to get the site terms:
    res_obj_path=path.expanduser(home+run_name+'/')+run_name+'_robj.pckl'
    
    #Open the station object:
    rfile=open(res_obj_path,'r')
    robj=pickle.load(rfile)
    rfile.close()
    
    #G##et out the stations and site terms:
    #Get unique station names:
    stations,u_stations_ind=unique(robj.sta,return_index=True)
    siteterms=around(array(robj.site_terms)[u_stations_ind],decimals=2)
    
    #Get standard deviation:
    siteterm_mean=mean(siteterms)
    siteterm_std=std(siteterms)
    
    #Write out stats to a file:
    sfile=open(statsfile,'w')
    sfile.write('Statistics for Run                  '+run_name+'\n')
    sfile.write('\n')
    sfile.write('Total Residual Mean:                '+str(mean_tot)+'\n')
    sfile.write('Total Residual Std Dev:             '+str(std_dev_tot)+'\n')
    sfile.write('\n')
    sfile.write('Event Residual Mean:                '+str(E_mean)+'\n')
    sfile.write('Event Residual Std Dev:             '+str(E_std_dev)+'\n')
    sfile.write('\n')
    sfile.write('Within-Event Residual Mean:         '+str(W_mean)+'\n')
    sfile.write('Within-Event Residual Std Dev:      '+str(W_std_dev)+'\n')
    sfile.write('\n')
    sfile.write('Site Term Residual Mean:            '+str(siteterm_mean)+'\n')
    sfile.write('Site Term Residual Std Dev:         '+str(siteterm_std)+'\n')
    sfile.write('\n')
    sfile.write('Path Term Residual Mean:            '+str(pterm_mean)+'\n')
    sfile.write('Path Term Residual Std Dev:         '+str(pterm_std)+'\n')
    sfile.write('\n')
    
    #Site terms:
    sfile.write('Site terms:\n')
    for site_i in range(len(siteterms)):
        sfile.write('\t'+stations[site_i]+':\t'+str(siteterms[site_i])+'\n')
    
    sfile.close()
    
    

#####
## Stuff for Mixed Effects ##
####
def plot_total(totalresid,home,run_name,resaxlim_mw):
    '''
    Plot total residuals given total residuals object from the mixed effects results
    Input:
        totalresid:      Object with total residuals including path
        home:            String with working home (i.e., /media/vsahakian/katmai/anza)
        run_name:        STring with run name (database/model combo)
        resaxlim_mw:      Array with r axis limits for the resid plot: [[mwmin,mwmax],[ymin,ymax]]
    Output:
    '''    

    import matplotlib.pyplot as plt
    from os import path
    
    #Get run directory:
    run_dir=path.expanduser(home+'/models/residuals/'+run_name+'/')
    
    totalfig = totalresid.plt_resids(run_name,resaxlim_mw)
    
    #figure names:
    totalfigname = run_dir + 'figs/' + run_name + '_total.png'
    totalfigpdf = run_dir + 'figs/pdfs/' + run_name + '_total.pdf'
    
    #Save:
    totalfig.savefig(totalfigname)
    totalfig.savefig(totalfigpdf)
    
    print 'Saved figures to %s and %s' % (totalfigname,totalfigpdf)
    

def plot_all(allresid,resaxlim_r,resaxlim_mw):
    '''
    Plot total residuals and path terms given  residuals object from the mixed effects results
    Input:
        allresid:        Object with all residuals including path
        resaxlim_r:      Array with axis limits for the resid plot: [[rmin,rmax],[ymin,ymax]]
        resaxlim_mw:     Array with axis limits for the resid plot: [[mwmin,mwmax],[ymin,ymax]]
    Output:
    '''
    
    import matplotlib.pyplot as plt
    
    ## First plot the total residuals:
    