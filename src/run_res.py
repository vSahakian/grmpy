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
        
        
def get_total_res(home,run_name,dbpath,modelpath,ffdf_flag,resaxlim):
    '''
    Get the total residuals for a database.
    Input:
        home:            STring with home path
        run_name:        String with name of the run
        dbpath:          String with path to a pickle database object
        modelpath:       String with path to a gmpe model pickle object
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
    from numpy import where
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
    d_predicted=gm.compute_model(model.m,model.rng,db.mw,db.r,db.ffdf,vs30,vref,mdep_ffdf)
    
    #Get residuals:
    total_residuals,mean_residual,std_dev=rcomp.total_residual(db,d_predicted)
    
    #Plot residuals and save plots:
    allresid=cdf.total_residuals(db.mw,total_residuals,mean_residual,std_dev)
    f1,f2=allresid.plt_resids(run_name,resaxlim)
    
    #Figure names:
    f1name=run_dir+'figs/'+run_name+'_total.png'
    f2name=run_dir+'figs/'+run_name+'_total_hist.png'
    #Save:
    f1.savefig(f1name)
    f2.savefig(f2name)
    
    return db.mw,allresid
    
    
    
    
def getEW_makeEvents(home,run_name,dbpath,modelpath,ffdf_flag):
    '''
    Get the Event and Within-Event Residuals, and store in Event objects
    Input:
        home:            STring with home path
        run_name:        String with name of the run
        dbpath:          String with path to a pickle database object
        modelpath:       String with path to a gmpe model pickle object
        ffdf_flag:       Flag for mag dependent ffdf.  0=off, 1=on
        resaxlim:        Array with axis limits for the resid plot: [[xmin,xmax],[ymin,ymax]]
    '''
    
    import numpy as np
    import cdefs as cdf
    import rescomp as rcomp
    from os import path
    import gmpe as gm
    import cPickle as pickle
    
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
    vref=760
    mdep_ffdf=ffdf_flag
    #Where are they 0?
    vs30_0ind=np.where(db.vs30==0)[0]
    #Get vs30 from database...
    vs30=db.vs30
    #Set 0 entries to vref:
    vs30[vs30_0ind]=vref
    
    #Now compute the predicted value of PGA...
    d_predicted=gm.compute_model(model.m,model.rng,db.mw,db.r,db.ffdf,vs30,vref,mdep_ffdf)

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
    E_std_dev=[]
    
    ##Zero out lists to use for indexing - use this resulting index, the station index,
    # to sort out the events for the station object later on...
    unique_stnums=np.unique(db.stnum)
    unique_stas=np.unique(db.sta)

    #Loop through the unique events, make each into an object, append to event list
    for i in range(len(unique_events)):
        unique_ind=np.where(unique_events[i]==db.evnum)[0]
        #Get the predicted data for this event only, for each recording
        #d_predicted_i is an array, with len = # of recordings for event):
        d_predicted_i=d_predicted[unique_ind]
    
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
        lat_i=db.lat[unique_ind]
        lon_i=db.lon[unique_ind]
        depth_i=db.depth[unique_ind]
        

        #The only one that is not directly from teh database is vs30; this is because
        #there were some 0's in it, so above that is compensated for by changing 
        #them to vref.
        vs30_i=vs30[unique_ind]
    
    
        #Make the event object:
        eventi=cdf.event(evnum_i,sta_i,stnum_i,ml_i,mw_i,pga_i,pgv_i,pga_pg_i,r_i,vs30_i,ffdf_i,md_ffdf_i,lat_i,lon_i,depth_i)
        
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
    
    return E_evnum,E_mw,E_residual,E_mean,E_std_dev
    
    
def sta_list(home,run_name,dbfile):
    '''
    '''
    
    import cdefs as cdf
    import numpy as np
    import cPickle as pickle
    import dread
    from os import path
    
    #Get directory for output:
    run_dir=path.expanduser(home+run_name+'/')
    so_dir=run_dir+'sta_objs/'
    fig_dir=run_dir+'figs/'
    
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
        lat=[]
        lon=[]
        depth=[]
        #Event residuals:
        E_residual=[]
        #Within-event residuals:
        W_residual=[]
        
        #Does this e
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
                lat.append(eventi.lat[sta_ev_ind])
                lon.append(eventi.lon[sta_ev_ind])
                depth.append(eventi.depth[sta_ev_ind])
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
        lat=np.array(lat)
        lon=np.array(lon)
        depth=np.array(depth)
        E_residual=np.array(E_residual)
        W_residual=np.array(W_residual)
        
        #Put into a station object...
        station_i=cdf.station(station_name_i[0],station_num_i,vs30[0],evnum,ml,mw,pga_pg,pga,pgv,ffdf,md_ffdf,lat,lon,depth,E_residual,W_residual)
        
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