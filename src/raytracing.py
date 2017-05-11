######Raytracing Module######
#VJS 8/2016

##Module to deal with raytracing stuff...
#Make the input files, plot results, etc.

def write_sourcein(home,run_name,veltype,lontype):
    '''
    Write a source.in file for fm3d raytracing
    VJS 8/2016
    Input:
        home:           String to home path for the main working directory
        run_name:       String to the name of the run
        veltype:        1=Vp, 2=Vs
        lontype:        Output lon type, assuming the input longitude in the database
                        is -118 etc. for west...
                        0=negative West (i.e., -118), 1=positive West (i.e., 241) 
    Output:
        source_veltype.in:     Writes out the source.in file for raytracing
    '''
    
    from os import path
    import cPickle as pickle
    import dread
    from numpy import str
    
    #Get the velocity type:
    if veltype==1:
        vm='Vp'
    elif veltype==2:
        vm='Vs'
    
    #Source is always local, for now:
    telesource=0
    #Number of paths always 1, since there's just one layer:
    npaths=1
    
    #Get the path to the run directory:
    run_dir=path.expanduser(home+run_name+'/')
    #Get the event object directory:
    eodir=run_dir+'event_objs/'
    #Get the path to the list of event objects:
    event_list_path=eodir+run_name+'.pckl'
    
    #Also get the path for the output sources.in file:
    raydir=run_dir+'raytracing/'
    sourcefile=raydir+'sources_'+vm+'.in'
    
    ##
    #Read in the list of event objects:
    eobjs=dread.read_obj_list(event_list_path)

    f=open(sourcefile,'w')
    
    #Get the number of sources:
    numsources=len(eobjs)
    #Get line to write:
    sline=str(numsources)+'\t\t\t\t number of sources\n'
    #Write to file:
    f.write(sline)
    
    #Loop over the event objects, and write each source (each event):
    for event_ind in range(len(eobjs)):
        #Get the event:
        eventi=eobjs[event_ind]
        
        if lontype==0:
            #Get the source info:
            #position
            sdepth=eventi.edepth[0]
            slat=eventi.elat[0]
            slon=eventi.elon[0]
            #number of paths - always =1 here, defined above
        elif lontype==1:
            sdepth=eventi.edepth[0]
            slat=eventi.elat[0]
            slon=eventi.elon[0]+360            
        
        #Number of sections on the path:  Always 1 here, since one layer
        nsections=1
        #Define path sections:
        psections=[0,2]
        #Velocity type:  #Specified as input (veltype)
        
        #Get lines to write:
        teleline=str(telesource)+'\t\t\t\t source is local/teleseismic (0/1)\n'
        posline=str(sdepth)+' '+str(slat)+' '+str(slon)+'\t\t\t position depth(km),lat(deg),lon(deg)\n'
        npathsline=str(npaths)+'\t\t\t\t number of paths from this source\n'
        nsectionsline=str(nsections)+'\t\t\t\t number of sections on the path\n'
        psectionsline=str(psections[0])+' '+str(psections[1])+'\t\t\t\t\t define the path sections\n'
        vtypeline=str(veltype)+'\t\t\t\t define the velocity type along the path\n'
        
        #Write them:
        f.write(teleline)
        f.write(posline)
        f.write(npathsline)
        f.write(nsectionsline)
        f.write(psectionsline)
        f.write(vtypeline)
        
    f.close()
        
    #Print statement
    print 'Wrote file '+sourcefile+' for velocity type '+vm+' to '+sourcefile


def write_receiverin(home,run_name,lontype):
    '''
    Write a receiver.in in file for fm3d raytracing
    Input:
        home:           String to the home path for the working directory
        run_name:       String for the run name (db/model combo)
        lontype:        Output lon type, assuming the input longitude in the database
                        is -118 etc. for west...
                        0=negative West (i.e., -118), 1=positive West (i.e., 241) 
    Output:
        receivers.in:   Receivers.in file for fm3d raytracing
    '''
    
    from os import path
    import cPickle as pickle
    import dread
    from numpy import str,ones,int64,shape
    
    #Get the path to the run directory:
    run_dir=path.expanduser(home+run_name+'/')
    #Get the event object directory:
    sodir=run_dir+'sta_objs/'
    #Get the path to the list of event objects:
    station_list_path=sodir+run_name+'.pckl'
    
    #Also get the path for the output sources.in file:
    raydir=run_dir+'raytracing/'
    receiverfile=raydir+'receivers.in'
    
    ##
    #Read in the list of station objects:
    sobjs=dread.read_obj_list(station_list_path)
    
    #Open the receiver file to write:
    f=open(receiverfile,'w')
    
    #Get the number of receivers:
    nrec=len(sobjs)
    #Get the line to write:
    nrecline=str(nrec)+'\t\t\t\t number of receivers\n'
    #Write line
    f.write(nrecline)
    
    #For each station, write the receiver information:
    for rec_ind in range(len(sobjs)):
        receiveri=sobjs[rec_ind]
        
        # Some of the residuals objects are made such that stlat is: array([[1,2,3,...,n]]).  If this is the case,
        #   shape(receiveri.stlat) = (1,n), and a [0] needs to be added on to indexing:
        if shape(receiveri.stlat)[0]==1:
            #Get info:
            #station info
            #Receiver elevation should be negative (station elevation as written is positive in database objects)
            rdepth=-1*receiveri.stelv[0][0]
            #get lat/lon - use lontype flag:
            if lontype==0:
                rlat=receiveri.stlat[0][0]
                rlon=receiveri.stlon[0][0]
            elif lontype==1:
                rlat=receiveri.stlat[0][0]
                rlon=receiveri.stlon[0][0]+360

            #number of paths to receiver
            npaths=len(receiveri.source_i[0])
            #The source of each path:
            sourcei=receiveri.source_i[0]
        
        
        # Otherwise, receiveri.stlat = array([1,2,3,...,n]), in which case
        #   shape(receiveri.stlat) = (n,), and don't need the extra 0...
        else:
            #Get info:
            #station info
            #Receiver elevation should be negative (station elevation as written is positive in database objects)
            rdepth=-1*receiveri.stelv[0]
            #get lat/lon - use lontype flag:
            if lontype==0:
                rlat=receiveri.stlat[0]
                rlon=receiveri.stlon[0]
            elif lontype==1:
                rlat=receiveri.stlat[0]
                rlon=receiveri.stlon[0]+360
            
            #number of paths to receiver
            npaths=len(receiveri.source_i)
            #The source of each path:
            sourcei=receiveri.source_i

        #The number of the source for each path - here it's always 1, since 
        #we're just raytracing through one layer...make them integers:
        npsource=ones((len(sourcei)))
        npsource=npsource.astype(int64)
        
        #Turn them into printable lines:
        rposline=str(rdepth)+' '+str(rlat)+' '+str(rlon)+'\t\t position depth(km),lat(deg),lon(deg)\n'
        nposline=str(npaths)+'\t\t\t\t number of paths to this receiver\n'
        
        #Turn the arrays into printable lines:
        #First the source number:
        for source in range(len(sourcei)):
            if source==0:
                sourceiline=str(sourcei[source])
            else:
                sourceiline=sourceiline+' '+str(sourcei[source])
        #end it:
        sourceiline=sourceiline+'\t\t the source of each path\n'
        
        #Now the number of the path in the source list:
        for pathi in range(len(npsource)):
            if pathi==0:
                npsourceline=str(npsource[pathi])
            else:
                npsourceline=npsourceline+' '+str(npsource[pathi])
        #end it:
        npsourceline=npsourceline+'\t\t\t number of the path in the source list\n'
        
        #Write them:
        f.write(rposline)
        f.write(nposline)
        f.write(sourceiline)
        f.write(npsourceline)
        
    #Close the file:
    f.close()
    
    #Print statememnt:
    print 'Wrote file '+receiverfile
    
        
    
    
###################################
#Function to read in a rayfile:

def parse_rayfile(rayfile):
    '''
    Parse rays.dat output file from fm3d
    Input:
        rayfile:        Path to the rays.dat file
    Output:
        path_list:      List of arrays, each array has path information
        receiver_id:    Array with id's of the receiver from receiver.in file
        source_id:      Array with id's of the source for that path from source.in file
        
    '''
    
    from numpy import array,r_,expand_dims
    
    #Load file
    f=open(rayfile,'r')
    
    #intialize things
    receiver_id=array([])
    source_id=array([])
    path_list=[]
    
    #Loop until end
    while True:
        
        #Read ray headers, this line has source and station IDs
        line=f.readline()
        
        #Have I recahed the end of the file?
        if line=='':
            break
            
        receiver_id=r_[receiver_id,int(line.split()[0])]
        source_id=r_[source_id,int(line.split()[1])]
        
        ##DEBUGGING...
        #print int(line.split()[0])
        #print int(line.split()[1])
        
        #Next line has number of elements in the path
        line=f.readline()
        Nelements=int(line.split()[0])
        
        #Initalize path coords variable
        path_coordinates=array([]) 
        
        #read path elements
        for kelement in range(Nelements):
            line=f.readline()
            element=array([float(line.split()[0]),float(line.split()[1]),float(line.split()[2])])
            #turn to row vector for appending
            element=expand_dims(element,0)
            if kelement==0:
                path_coordinates=element
            else:
                path_coordinates=r_[path_coordinates,element]
            
        #Append to list
        path_list.append(path_coordinates)

    f.close()
    
    print 'read rayfile'

    return path_list,receiver_id,source_id
        
        

#######
#Function to extract rayfile info, sort, convert to depth/lat/lon, and store
#in a residuals object
def store_rayinfo(rfile_in,rfile_out,rayfile,veltype,lon_conversion):
    '''
    Extract rayfile info, sort, convert to depth/lat/lon, and store into
    the residuals object
    Input:
        rfile_in:       String with input residual object location
        frile_out:      String with output residual object location
        rayfile:        String with the path to the rayfile
        veltype:        Velocity type Vp=1, Vs=2
        lon_conversion: Flag to convert longitude.  0=none; 1=+West to -West; 2=-West to +West
    Output:
        residuals object
    '''
    
    from os import path
    import cPickle as pickle
    from numpy import where,append,degrees,array
    
    print rfile_in
    
    #Load the residuals object:
    rfile=open(rfile_in,'r')
    robj=pickle.load(rfile)
    rfile.close()
    
    #Get the raytracing info:
    path_list,rec_id,src_id=parse_rayfile(rayfile)
    
    #Zero out the sorted path list:
    ray_depth=[]
    ray_lat=[]
    ray_lon=[]
    
    #Find which raypath corresponds to each recording in the database/residuals object:
    for recording_i in range(len(robj.source_i)):
        source_i=robj.source_i[recording_i]
        receiver_i=robj.receiver_i[recording_i]
        
        print robj.source_i[recording_i]
        
        #Find which raypath entry corresponds to this recording:
        #print source_i
        #print src_id
        #print receiver_i
        #print rec_id
        ray_ind=where((source_i==src_id) & (receiver_i==rec_id))[0]
        print ray_ind
        
        #Convert the ray info to depth/lat degrees/lon degrees:
        ray_depth_i=path_list[ray_ind][:,0]-6371   #minus the radius of the earth
        ray_lat_i=degrees(path_list[ray_ind][:,1])
        #convert longitude?
        if lon_conversion==0:
            ray_lon_i=degrees(path_list[ray_ind][:,2])
        elif lon_conversion==1:
            ray_lon_i=degrees(path_list[ray_ind][:,2])-360
        elif lon_conversion==2:
            ray_lon_i=degrees(path_list[ray_ind][:,2])+360
        
        #Add this info to the sorted path list:
        ray_depth.append(ray_depth_i)
        ray_lat.append(ray_lat_i)
        ray_lon.append(ray_lon_i)

    
    #Store back into a residuals object
    if veltype==1:
        #Then it's vp.  Write to the Vp part of the db.
        robj.add_vp_paths(ray_depth,ray_lat,ray_lon)
    elif veltype==2:
        #Then it's vs.  Write the Vs part of the db.
        robj.add_vs_paths(ray_depth,ray_lat,ray_lon)
        
        
    #Save the file again...
    rfile=open(rfile_out,'w')
    pickle.dump(robj,rfile)
    rfile.close()

    
##############    
def plot_rays(home,run_name,dbname,veltype,view,axlims,stations,events,by_path,mymap,faultfile,residfilepath='none'):
    '''
    Plot the raypaths and save the png and pdf figures
    VJS 8/2016
    Input:
        home:           String with the home directory for the project
        run_name:       String with the run name combo of db and model
        dbname:         String with the basename of the residuals object
        veltype:        Vp/Vs (1/2)
        view:           Plot view; map/lat vs depth/lon vs depth, 0/1/2
        axlims:         Axis limits [[xmin, xmax],[ymin,ymax]]
        stations:       Plot stations?  no/yes = 0/1
        events:         Plot events?  no/yes = 0/1
        by_path:        Color raypaths by path?  no/yes=0/1
        mymap:          String with the python colormap to use
        faultfile:      String to the path of the pckl file with faults to plot
        residfilepath:  String with path to the residuals file, if it doesn't follow the typical _raydat.pckl format
    Output: 
        figure          Prints a png and pdf version of the figure to the run fig directory
    '''
    
    from os import path
    import cPickle as pickle
    
    #Get the run directory:
    run_dir=path.expanduser(home+run_name+'/')
    
    #And the figure directories:
    fig_dir=run_dir+'figs/'
    pdf_dir=fig_dir+'pdfs/'
    
    #Get the residuals object:
    if residfilepath=='none':
        residfile_old=run_dir+dbname+'_robj.pckl'
        #Get the new name to add on raydat:
        rbase=residfile_old.split('.pckl')
        residfile=rbase[0]+'_raydat.pckl'
    else:
        residfile=residfilepath
    
    #Load the residuals object:
    rfile=open(residfile,'r')
    robj=pickle.load(rfile)
    rfile.close()
    
    figure=robj.plot_raypaths(veltype,view,axlims,stations,events,by_path,mymap,faultfile)
    
    #Save the figures:
    #Get the figure name:
    pngfile=fig_dir+run_name+'_vtype'+str(veltype)+'_view'+str(view)+'_sta'+str(stations)+'_ev'+str(events)+'.png'
    pdffile=pdf_dir+run_name+'_vtype'+str(veltype)+'_view'+str(view)+'_sta'+str(stations)+'_ev'+str(events)+'.pdf'
    #Save png:
    figure.savefig(pngfile)
    #save pdf:
    figure.savefig(pdffile)
    
    
    
##############    
def plot_rays_together(home,run_name,dbname,veltype,view,axlims,stations,events,by_path,mymap,faultfile,residfilepath='none'):
    '''
    Plot the raypaths and save the png and pdf figures
    VJS 8/2016
    Input:
        home:           String with the home directory for the project
        run_name:       String with the run name combo of db and model
        dbname:         String with the basename of the residuals object
        veltype:        Array or list with the velocity types to plot: 0 = Vp, 1 = Vs
        view:           Array or list with view types to plot: map=0, lat vs depth=1, lon vs depth=2
        axlims:         Axis limits for each view type, 1 - n: [[[xmin1, xmax1],[ymin1,ymax1]],...,[[xmin_n, xmax_n],[ymin_n,ymax_n]]]
        stations:       Plot stations?  no/yes = 0/1
        events:         Plot events?  no/yes = 0/1
        by_path:        Color raypaths by path?  no/yes=0/1
        mymap:          String with the python colormap to use
        faultfile:      String to the path of the pckl file with faults to plot
        residfilepath:  String with path to the residuals file, if it doesn't follow the typical _raydat.pckl format
    Output: 
        figure          Prints a png and pdf version of the figure to the run fig directory
    '''
    
    from os import path
    import cPickle as pickle
    from numpy import str
    
    #Get the run directory:
    run_dir=path.expanduser(home+run_name+'/')
    
    #And the figure directories:
    fig_dir=run_dir+'figs/'
    pdf_dir=fig_dir+'pdfs/'
    
    #Get the residuals object:
    if residfilepath=='none':
        residfile_old=run_dir+dbname+'_robj.pckl'
        #Get the new name to add on raydat:
        rbase=residfile_old.split('.pckl')
        residfile=rbase[0]+'_raydat.pckl'
    else:
        residfile=residfilepath
    
    #Load the residuals object:
    rfile=open(residfile,'r')
    robj=pickle.load(rfile)
    rfile.close()
    
    for vtype_i in range(len(veltype)):
        for view_i in range(len(view)):
                figure=robj.plot_raypaths(veltype[vtype_i],view[view_i],axlims[view_i],stations,events,by_path,mymap,faultfile)
    
                #Save the figures:
                #Get the figure name:
                pngfile=fig_dir+run_name+'_vtype'+str(veltype[vtype_i])+'_view'+str(view[view_i])+'_sta'+str(stations)+'_ev'+str(events)+'.png'
                pdffile=pdf_dir+run_name+'_vtype'+str(veltype[vtype_i])+'_view'+str(view[view_i])+'_sta'+str(stations)+'_ev'+str(events)+'.pdf'
                #Save png:
                figure.savefig(pngfile)
                #save pdf:
                figure.savefig(pdffile)
                
                printstring = 'plotted and saved raypaths for view ' + str(view[view_i]) + ', and velocity type ' + veltype[vtype_i]
                print printstring
    
    
    
    
######
def plot_rays_cutoffval(home,run_name,dbname,veltype,view,axlims,stations,events,mymap,faultfile,cutoff_val,residfilepath='none'):
    '''
    Plot the raypaths and save the png and pdf figures
    VJS 8/2016
    Input:
        home:           String with the home directory for the project
        run_name:       String with the run name combo of db and model
        dbname:         STring with the basename of the residuals object
        veltype:        Vp/Vs (1/2)
        view:           Plot view; map/lat vs depth/lon vs depth, 0/1/2
        axlims:         Axis limits [[xmin, xmax],[ymin,ymax]]
        stations:       Plot stations?  no/yes = 0/1
        events:         Plot events?  no/yes = 0/1
        by_path:        Color raypaths by path?  no/yes=0/1
        mymap:          String with the python colormap to use
        faultfile:      String to the path of the pckl file with faults to plot
        cutoff_val:     Value above/below which to color raypath; otherwise the path is gray (plots if abs(path_term)>=cutoff_val)
        residfilepath:  String with path to the residuals file, if it doesn't follow the typical _raydat.pckl format

    Output: 
        figure          Prints a png and pdf version of the figure to the run fig directory
    '''
    
    from os import path
    import cPickle as pickle
    
    #Get the run directory:
    run_dir=path.expanduser(home+run_name+'/')
    
    #And the figure directories:
    fig_dir=run_dir+'figs/'
    pdf_dir=fig_dir+'pdfs/'
    
    
    #Get the residuals object:
    if residfilepath=='none':
        residfile_old=run_dir+dbname+'_robj.pckl'
        #Get the new name to add on raydat:
        rbase=residfile_old.split('.pckl')
        residfile=rbase[0]+'_raydat.pckl'
    else:
        residfile=residfilepath
    
    #Load the residuals object:
    rfile=open(residfile,'r')
    robj=pickle.load(rfile)
    rfile.close()
    
    figure=robj.plot_raypaths_cutoffval(veltype,view,axlims,stations,events,mymap,faultfile,cutoff_val)
    
    #Save the figures:
    #Get the figure name:
    pngfile=fig_dir+run_name+'_vtype'+str(veltype)+'_view'+str(view)+'_sta'+str(stations)+'_ev'+str(events)+'_cval'+str(cutoff_val)+'.png'
    pdffile=pdf_dir+run_name+'_vtype'+str(veltype)+'_view'+str(view)+'_sta'+str(stations)+'_ev'+str(events)+'_cval'+str(cutoff_val)+'.pdf'
    #Save png:
    figure.savefig(pngfile)
    #save pdf:
    figure.savefig(pdffile)
    
    
#######
def plot_3d_raypaths(home,run_name,dbname,vtype,stations,events,axlims,colormap,faultfile):
    '''
    Plot the 3d raypaths for a given object
    Input:
        home:           String with the path to the project home
        run_name:       String with the database/model combination
        dbname:         STring with the basename of the residuals object
        vtype:          Velocity type to plot: 1/2 = Vp/Vs
        stations:       Plot stations?  0/1=no/yes
        events:         Plot events?    0/1=no/yes
        axlims:         Axis limits: [[xmin,xmax],[ymin,ymax],[zmin,zmax]]
        colormap:       STring with colormap to plot, i.e. 'jet'
        faultfile:      String with path to the faults stored in a pckl file
    Output:
        figure3d:       Figure in 3d
    '''
    
    from os import path
    import cPickle as pickle
    import matplotlib.pyplot as plt
    
    #Get the run directory:
    run_dir=path.expanduser(home+run_name+'/')
    
    #Residuals file with raypaths:
    residpath=run_dir+dbname+'_robj_raydat.pckl'

    #Open the object:
    rfile=open(residpath,'r')
    robj=pickle.load(rfile)
    rfile.close()
    
    #plot:
    figure3d=robj.plot_raypaths_3d(vtype,stations,events,axlims,colormap,faultfile)
    
    #return:
    plt.show(figure3d)
    
    return figure3d
    
###########
def pull_rays_frombox_forstation(r_object,bounds,station_name,ray_type,mymap):
    '''
    Input:
        r_object:           Residuals object with raypaths stored
        bounds:             Bounds for box of events: [w,e,n,s]
        station_name:       String with station name
        ray_type:           Type of ray: 'Vp' or 'Vs'
        mymap:              Colormap string
    Output:
    '''
    
    import cPickle as pickle
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import path
    import matplotlib.colors as colors
    import matplotlib.cm as cm
    import matplotlib.patches as patches

    # What raytype is being used?
    if ray_type == 'Vp':
        ray_lon = r_object.vp_lon
        ray_lat = r_object.vp_lat
        ray_depth = r_object.vp_depth
    elif ray_type == 'Vs':
        ray_lon = r_object.vs_lon
        ray_lat = r_object.vs_lat
        ray_depth = r_object.vs_depth
        
    # Define boundaries for events to keep:
    w = bounds[0]
    e = bounds[1]
    n = bounds[2]
    s = bounds[3]
        
    
    # The station to select is the name given:
    select_sta = station_name
    
    # Find the indices of this station within the residuals object:
    select_ind = np.where(r_object.sta==select_sta)[0]
    
    # Get station, event, and path info for that station:
    select_stlon = r_object.stlon[select_ind][0]
    select_stlat = r_object.stlat[select_ind][0]
    select_stelv = r_object.stelv[select_ind][0]
    
    # Events recorded:
    select_evnum = r_object.evnum[select_ind]
    select_elon = r_object.elon[select_ind]
    select_elat = r_object.elat[select_ind]
    select_edepth = r_object.edepth[select_ind]
    
    # raypaths:
    select_ray_depth = np.array(ray_depth)[select_ind]
    select_ray_lon = np.array(ray_lon)[select_ind]
    select_ray_lat = np.array(ray_lat)[select_ind]
    
    # path terms:
    select_path_terms = r_object.path_terms[select_ind]
        
    #Get mean and std of dataset for plotting ray color:
    mean_pterm=np.mean(select_path_terms)
    std_pterm= np.std(select_path_terms)
    #Set hte colorscale to cover 97% of the data:
    cmin=-3*std_pterm
    cmax=3*std_pterm   

    ######################
    ##Plot:
    #Get colormap
    #Make colormap:
    colormap_pterm=plt.get_cmap(mymap)
    #Make a normalized colorscale
    cNorm=colors.Normalize(vmin=cmin, vmax=cmax)
    #Apply normalization to colormap:
    scalarMap=cm.ScalarMappable(norm=cNorm, cmap=colormap_pterm)
    
    #Make a fake contour plot for the colorbar:
    dummyplot = plt.figure()
    Z=[[0,0],[0,0]]
    levels=np.arange(cmin,cmax,0.01)
    c=plt.contourf(Z, levels, cmap=colormap_pterm)
    plt.close(dummyplot)

    # Plot all rays recorded at this station:
    all_rays = plt.figure()   
    ax = all_rays.add_subplot(111, aspect='equal')
    for ray in range(len(select_ind)):
        
        #Assign color to path term:
        colorVal = scalarMap.to_rgba(select_path_terms[ray])
        
        # Plot:
        ax.plot(select_ray_lon[ray],select_ray_lat[ray],color=colorVal)  
        
    ax.add_patch(
        patches.Rectangle(
            (w, s),
            (e - w),
            (n - s),
            fill=False      # remove background
        )
    ) 
    
    # Set colorbar
    cb = plt.colorbar(c)
    cb.set_label('Path term (ln residual)')
    # Set title:
    maptitle = 'Map view plot of station %s, mean path term is %.2f' % (select_sta,np.mean(select_path_terms))
    plt.title(maptitle)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    
    
    ### Now extract data from the box:    
    # Make the path:
    keeppath = path.Path([[w,s],[e,s],[e,n],[w,n],[w,s]])
    
    # INitiate array for keeping the events within the path:
    keep_event_ind=np.array([]).astype('int')
    
    # Get the rays at this station within the box
    for rayi in range(len(select_ind)):
        if keeppath.contains_point([select_elon[rayi],select_elat[rayi]]):
            keep_event_ind = np.r_[keep_event_ind,rayi]
    
    
    # Plot now only the rays fromt hat box, recorded on this station
    #Make a fake contour plot for the colorbar:
    dummyplot = plt.figure()
    Z=[[0,0],[0,0]]
    levels=np.arange(cmin,cmax,0.01)
    c=plt.contourf(Z, levels, cmap=colormap_pterm)
    plt.close(dummyplot)    
    
    box_rays_map = plt.figure()
    for ray in range(len(keep_event_ind)):
        colorVal = scalarMap.to_rgba(select_path_terms[keep_event_ind][ray])
        plt.plot(select_ray_lon[keep_event_ind[ray]],select_ray_lat[keep_event_ind[ray]],color=colorVal)
    
    # Set colorbar
    cb = plt.colorbar(c)
    cb.set_label('Path term (ln residual)')
    # Set title:
    maptitle = 'Map view plot of box, mean path term is %.2f' % np.mean(select_path_terms[keep_event_ind])
    plt.title(maptitle)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
       
       
    ######   
    # Plot them now in cross section with longitude   
    #Make a fake contour plot for the colorbar:
    dummyplot = plt.figure()
    Z=[[0,0],[0,0]]
    levels=np.arange(cmin,cmax,0.01)
    c=plt.contourf(Z, levels, cmap=colormap_pterm)
    plt.close(dummyplot)
    
    box_rays_londep = plt.figure()
    for ray in range(len(keep_event_ind)):
        colorVal = scalarMap.to_rgba(select_path_terms[keep_event_ind][ray])
        plt.plot(select_ray_lon[keep_event_ind[ray]],select_ray_depth[keep_event_ind[ray]],color=colorVal)
        
    # Set colorbar
    cb = plt.colorbar(c)
    cb.set_label('Path term (ln residual)')    
    # Set title:
    lontitle = 'Longitude plot of box, mean path term is %.2f' % np.mean(select_path_terms[keep_event_ind])
    plt.title(lontitle)   
    plt.xlabel('Longitude')
    plt.ylabel('Depth (km)') 
                
                                
    # Return the figures, and the indices to keep - station indices, and keep event indices together:
    return select_ind, keep_event_ind, all_rays, box_rays_map, box_rays_londep
    
    
def save_box_data(text_dir,fig_dir,pdf_dir,bounds,select_sta,bound_name,all_rays_fig,box_rays_map,box_rays_londep):
    '''
    Input: 
        text_dir:                   String with text directory for bounds info
        fig_dir:                    String with png directory for figures
        pdf_dir:                    String with pdf directory for figures
        boudns:                     List:  [w,e,n,s]
        select_sta:                 String with station name
        bound_name:                 String with box bound name, i.e. '3a'
        all_rays_fig:               PYthon figure with map view of all rays per station
        box_rays_map:               Python figure with map view of all rays at station in box
        box_rays_londep:            Python figure with lon/depth view of all rays at station in box
    Output:
        Saves figures to png, pdf, and txt files.
    '''
    import matplotlib.pyplot as plt
    
    ###  Save figures, save text file with bounds:
    ## text file
    combo_name = select_sta + '_'+ bound_name
    header = ['w','e','n','s']
    
    ## bounds text file
    txtfile = text_dir + combo_name + '.txt'
    f = open(txtfile,'w')
    for corneri in range(len(header)):
        writeline = '%s \t %.5f \n' % (header[corneri],bounds[corneri])
        f.write(writeline)
    f.close()
    
    ## all rays at station:
    png_station = fig_dir + combo_name + '_sta.png'
    pdf_station = pdf_dir + combo_name + '_sta.pdf'
    
    all_rays_fig.savefig(png_station)
    all_rays_fig.savefig(pdf_station)

    ## rays in box, map view
    png_map = fig_dir + combo_name + '_boundsmap.png'
    pdf_map = pdf_dir + combo_name + '_boundsmap.pdf'
    
    box_rays_map.savefig(png_map)
    box_rays_map.savefig(pdf_map)

    ## rays in box, lon view
    png_cross = fig_dir + combo_name + '_boundscross.png'
    pdf_cross = pdf_dir + combo_name + '_boundscross.pdf'
    
    box_rays_londep.savefig(png_cross)
    box_rays_londep.savefig(pdf_cross)