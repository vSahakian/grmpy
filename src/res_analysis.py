######Residual Analysis######
##VJS 6/2016

##Interpolate rays through a material model##
def interpolate_rays(residobj,materialobj,interptype,rayflag):
    '''
    Interpolate rays through a material model
    Input:
        residobj:           Object with residuals and ray position information
        materialobj:        Object with material model.  Depth should be positive for input.
        interptype:         String with type of interpolation from rbf, i.e., 'linear'
        rayflag:            Flag for type of ray to interpolate.  0=Vp/1=Vs
    Output:
        ray_data:           List of arrays, each of same length as corresponding ray.
                            Arrays contain values of model at those raypath positions.
    '''
    
    from numpy import meshgrid,zeros,where
    from scipy.interpolate import Rbf

    #Get the actual grid x and y values:
    grid_x=materialobj.x
    grid_y=materialobj.y
    grid_z=-materialobj.z

    #Turn these into a grid, prepping for ravel for interpolation:
    gX,gY,gZ=meshgrid(grid_x,grid_y,grid_z,indexing='ij')

    #Make column vectors to put into the interpolation:
    columnx=gX.ravel()
    columny=gY.ravel()
    columnz=gZ.ravel()
    #Data to interpolate - transpose so the x is first, etc.:
    data=materialobj.materials.transpose(2,1,0).ravel()  
    
    #Make interpolator
    print 'Making interpolator...'
    interpolator = Rbf(columnx, columny, columnz, data,function=interptype)

    #Get the values to interpolate., then interpolate
    #First, which value are you doing this for?
    if rayflag==0:
        #Initialize ray_data array:
        ray_data=[]
        
        print 'Interpolating P rays'
        #Run in a loop for each ray:
        for ray_i in range(len(residobj.vp_lon)):
            ray_x_i=residobj.vp_lon[ray_i]
            ray_y_i=residobj.vp_lat[ray_i]
            ray_z_i=residobj.vp_depth[ray_i]
                
            #Interpolate:
            vp_i=interpolator(ray_x_i,ray_y_i,ray_z_i)
    
            #Append to the list:
            ray_data.append(vp_i)
            
            #Print position:
            if ray_i % 1000 ==0: print ray_i
            
            
    elif rayflag==1:
        #Initialize ray_data array:
        ray_data=[]
        
        print 'Interpolating S rays'
        for ray_i in range(len(residobj.vs_lon)):
            ray_x_i=residobj.vs_lon[ray_i]
            ray_y_i=residobj.vs_lat[ray_i]
            ray_z_i=residobj.vs_depth[ray_i]
                
            #Interpolate:
            vs_i=interpolator(ray_x_i,ray_y_i,ray_z_i)
    
            #Append to the list:
            ray_data.append(vs_i)            
            
            #Print position:
            if ray_i % 1000 ==0: print ray_i
            
    #Return info:
    return ray_data
        
        
############
##Compute indices##
#Path integral:
def compute_pathintegral(ray_vals,materialobject,normalize_flag):
    '''
    Compute a path integral index through a material model
    Input:
        ray_vals:               List of arrays with the values of the model for each ray
        materialobject:         Object with material model
        normalize_flag:         Normalize path integral? 0=no/1=yes
    Output:
        pathintegral_index:     Array with the path integral index
    '''
    
    from numpy import max,sum,array
    
    #Find the maximum value in the materials object to use for normalization integral
    if normalize_flag==1:
        max_material_val=max(materialobject.materials)
        
    #Initialize the index list:
    pathintegral_index=[]
    
    #Loop over rays:
    for ray_i in range(len(ray_vals)):
        #Get the values for this ray:
        ray_val_i=ray_vals[ray_i]
        #Get the length of this ray:
        ray_i_len=len(ray_val_i)
        
        if normalize_flag==1:
            #Get the "path integral" of the maximum value of hte array by
            #   multiplying hte maximum value by the number of points on this array:
            max_val_integral=ray_i_len*max_material_val
        
        #Sum the values of this ray:
        path_integral_i=sum(ray_val_i)
        
        #Set the value of the index for this ray:
        if normalize_flag==0:
            path_integral_index_i=path_integral_i
        elif normalize_flag==1:
            path_integral_index_i=path_integral_i/max_val_integral
            
        #Append to the index list:
        pathintegral_index.append(path_integral_index_i)
    
    #After appending all, turn it into an array:
    pathintegral_index=array(pathintegral_index)
        
    #Return the array:
    return pathintegral_index
    
#####
#Path integral:
def compute_devpathintegral(ray_vals,materialobject,normalize_flag):
    '''
    Compute a gradient of the path integral index through a material model
    Input:
        ray_vals:               List of arrays with the values of the model for each ray
        materialobject:         Object with material model
        normalize_flag:         Normalize path integral? 0=no/1=yes
    Output:
        dpathintegral_index:     Array with the gradient path integral index
    '''
    
    from numpy import max,sum,array,gradient
    
    ##Find the maximum value in the materials object to use for normalization integral
    #if normalize_flag==1:
    #    max_material_val=max(materialobject.materials)
        
    #Initialize the index list:
    dpathintegral_index=[]
    
    #Loop over rays:
    for ray_i in range(len(ray_vals)):
        #Get the values for this ray:
        ray_val_i=ray_vals[ray_i]

        #Get the length of this ray:
        ray_i_len=len(ray_val_i)
        
        #if normalize_flag==1:
        #    #Get the "path integral" of the maximum value of hte array by
        #    #   multiplying hte maximum value by the number of points on this array:
        #    max_val_integral=ray_i_len*max_material_val
        
        #Get the gradient along the ray:
        gradient_i=gradient(ray_val_i)
        
        #Sum the values of this ray:
        dpath_integral_i=sum(abs(gradient_i))
        
        #Set the value of the index for this ray:
        if normalize_flag==0:
            dpath_integral_index_i=dpath_integral_i
        #elif normalize_flag==1:
        #    path_integral_index_i=path_integral_i/max_val_integral
            
        #Append to the index list:
        dpathintegral_index.append(dpath_integral_index_i)
    
    #After appending all, turn it into an array:
    dpathintegral_index=array(dpathintegral_index)
        
    #Return the array:
    return dpathintegral_index
            
            
            
 #####################
 #####################
 ##Plots and Correlations##
 
def plot_pterms(home,run_name,robj,index,axlims):
    '''
    Plot path terms vs. some index, with path terms on the x-axis.
    Input:
        home:           String with path to the home directory for project
        run_name:       String with database/inversion combo for residuals
        robj:           Residuals object with indices
        index:          STring with index to plot: 'ind'_raytype_modeltype_indextype
                            raytype=p/s; modeltype=vp/vs/vpvs/qp/qs/qpqs;
                            indextype=pathint,normpathint,gradpathint
                            i.e., 'ind_p_vs_normpathint'
        axlims:         Axis limits: [[xmin,xmax],[ymin,ymax]]
    Output:
        Saves plots to home/run_name/figs/pathterm_index
    '''
    from os import path
    import matplotlib.pyplot as plt
    from scipy.stats.stats import pearsonr
    from numpy import array
    
    #Get run directory, and figure directory:
    run_dir=path.expanduser(home+run_name+'/')
    fig_dir=run_dir+'figs/'
    pdf_dir=fig_dir+'pdfs/'
    
    #Get a string for the plot name and title:
    figbasename='pathterm_'+index
    #Index type...
    indtype=index.split('_')[-1]
    if indtype=='pathint':
        indname='Path Integral'
    elif indtype=='normpathint':
        indname='Normalized Path Integral'
    elif indtype=='gradpathint':
        indname='Path Integral of the gradient'
        
    #What to plot?
    x=robj.path_terms
    y=getattr(robj,index)
        
    #get correlation coefficient:
    pcoeff,tails=pearsonr(x,y)
    pcoeff=round(pcoeff,2)
    
    #Title:
    ptitle='Plot of path terms vs. '+indname+'\n Pearson coefficient: '+str(pcoeff)
    
    #Initiate figure:
    f1=plt.figure()
    
    #Plot it:
    color=array([102,139,139])/255.0
    plt.scatter(x,y,edgecolors=color,facecolors='none',lw=0.9)
    
    #Set axis limits and labels:
    plt.xlim(axlims[0])
    plt.ylim(axlims[1])
    
    plt.xlabel('Path term (ln residual)')
    plt.ylabel(indname)
    plt.title(ptitle)
    
    #Save figure:
    pngname=fig_dir+index+'.png'
    pdfname=pdf_dir+index+'.pdf'
    plt.savefig(pngname)
    plt.savefig(pdfname)
          
    #Return
    return f1
    
##
def plot_pathterms_colored(home,run_name,robj,index,axlims,color_by,cvals,mymap):
    '''
    Plot path terms vs. some index, with path terms on the x-axis.
    Input:
        home:           String with path to the home directory for project
        run_name:       String with database/inversion combo for residuals
        robj:           Residuals object with indices
        index:          STring with index to plot: 'ind'_raytype_modeltype_indextype
                            raytype=p/s; modeltype=vp/vs/vpvs/qp/qs/qpqs;
                            indextype=pathint,normpathint,gradpathint
                            i.e., 'ind_p_vs_normpathint'
        axlims:         Axis limits: [[xmin,xmax],[ymin,ymax]]
        color_by:       String with variable to color by: 'r', or 'mw'
        cvals:          Colorbar limits: [cmin,cmax]
        mymap:          String with colormap to use (i.e., 'jet')
    Output:
        Saves plots to home/run_name/figs/pathterm_index
    '''
    from os import path
    import matplotlib.pyplot as plt
    from scipy.stats.stats import pearsonr
    from numpy import array,arange
    import matplotlib.colors as colors
    import matplotlib.cm as cm
    
    #Get run directory, and figure directory:
    run_dir=path.expanduser(home+run_name+'/')
    fig_dir=run_dir+'figs/'
    pdf_dir=fig_dir+'pdfs/'
    
    #Get a string for the plot name and title:
    figbasename='pathterm_'+index+'_'+color_by
    #Index type...
    indtype=index.split('_')[-1]
    if indtype=='pathint':
        indname='Path Integral'
    elif indtype=='normpathint':
        indname='Normalized Path Integral'
    elif indtype=='gradpathint':
        indname='Path Integral of the gradient'
        
    #What to plot?
    x=robj.path_terms
    y=getattr(robj,index)
        
    #get correlation coefficient:
    pcoeff,tails=pearsonr(x,y)
    pcoeff=round(pcoeff,2)
    
    #Title:
    ptitle='Plot of path terms vs. '+indname+'\n Pearson coefficient: '+str(pcoeff)
    
    #Get colormap
    #Make colormap:
    colormap_mwdist=plt.get_cmap(mymap)
    #Make a normalized colorscale
    cmin=cvals[0]
    cmax=cvals[1]
    cNorm=colors.Normalize(vmin=cmin, vmax=cmax)
    #Apply normalization to colormap:
    scalarMap=cm.ScalarMappable(norm=cNorm, cmap=colormap_mwdist)
    
    #Make a fake contour plot for the colorbar:
    Z=[[0,0],[0,0]]
    levels=arange(cmin,cmax,0.01)
    c=plt.contourf(Z, levels, cmap=colormap_mwdist)   
    
    #Assign values to colormap
    colorVal = scalarMap.to_rgba(getattr(robj,color_by))
    
    #Plot:
    f1=plt.figure()
    plt.scatter(x,y,facecolors='none',edgecolors=colorVal,lw=0.5)
    
    #Add colorbar:
    cb=plt.colorbar(c)
    if color_by=='r':
        cbarlabel='Distance (km)'
    elif color_by=='mw':
        cbarlabel='M'
    cb.set_label(cbarlabel)
    
    #Set axis limits and labels:
    plt.xlim(axlims[0])
    plt.ylim(axlims[1])
    
    plt.xlabel('Path term (ln residual)')
    plt.ylabel(indname)
    plt.title(ptitle)
    
    #Save figure:
    pngname=fig_dir+index+'_'+color_by+'.png'
    pdfname=pdf_dir+index+'_'+color_by+'.pdf'
    plt.savefig(pngname)
    plt.savefig(pdfname)
    
    return f1
    
##
#Other terms...
##
def plot_terms_colored(home,run_name,robj,term,index,axlims,color_by,cvals,mymap):
    '''
    Plot path terms vs. some index, with path terms on the x-axis.
    Input:
        home:           String with path to the home directory for project
        run_name:       String with database/inversion combo for residuals
        robj:           Residuals object with indices
        term:           String with term to plot: 'site_terms','W_residual','E_residual'
        index:          STring with index to plot: 'ind'_raytype_modeltype_indextype
                            raytype=p/s; modeltype=vp/vs/vpvs/qp/qs/qpqs;
                            indextype=pathint,normpathint,gradpathint
                            i.e., 'ind_p_vs_normpathint'
        axlims:         Axis limits: [[xmin,xmax],[ymin,ymax]]
        color_by:       String with variable to color by: 'r', or 'mw'
        cvals:          Colorbar limits: [cmin,cmax]
        mymap:          String with colormap to use (i.e., 'jet')
    Output:
        Saves plots to home/run_name/figs/pathterm_index
    '''
    from os import path
    import matplotlib.pyplot as plt
    from scipy.stats.stats import pearsonr
    from numpy import array,arange
    import matplotlib.colors as colors
    import matplotlib.cm as cm
    
    #Get run directory, and figure directory:
    run_dir=path.expanduser(home+run_name+'/')
    fig_dir=run_dir+'figs/'
    pdf_dir=fig_dir+'pdfs/'
    
    #Get a string for the plot name and title:
    figbasename=term+index+'_'+color_by
    #Index type...
    indtype=index.split('_')[-1]
    if indtype=='pathint':
        indname='Path Integral'
    elif indtype=='normpathint':
        indname='Normalized Path Integral'
    elif indtype=='gradpathint':
        indname='Path Integral of the gradient'
        
    #Term type:
    if term=='site_terms':
        termname='Site Term (ln residual)'
        termtitle='Site Term'
    elif term=='W_residual':
        termname='Within Event residual (ln residual)'
        termtitle='Within-Event Residual'
    elif term=='E_residual':
        termname='Event residual (ln residual)'
        termtitle='Event Residual'
        
    #What to plot?
    x=getattr(robj,term)
    y=getattr(robj,index)
        
    #get correlation coefficient:
    pcoeff,tails=pearsonr(x,y)
    pcoeff=round(pcoeff,2)
    
    #Title:
    ptitle='Plot of '+termtitle+' vs. '+indname+'\n Pearson coefficient: '+str(pcoeff)
    
    #Get colormap
    #Make colormap:
    colormap_mwdist=plt.get_cmap(mymap)
    #Make a normalized colorscale
    cmin=cvals[0]
    cmax=cvals[1]
    cNorm=colors.Normalize(vmin=cmin, vmax=cmax)
    #Apply normalization to colormap:
    scalarMap=cm.ScalarMappable(norm=cNorm, cmap=colormap_mwdist)
    
    #Make a fake contour plot for the colorbar:
    Z=[[0,0],[0,0]]
    levels=arange(cmin,cmax,0.01)
    c=plt.contourf(Z, levels, cmap=colormap_mwdist)   
    
    #Assign values to colormap
    colorVal = scalarMap.to_rgba(getattr(robj,color_by))
    
    #Plot:
    f1=plt.figure()
    plt.scatter(x,y,facecolors='none',edgecolors=colorVal,lw=0.5)
    
    #Add colorbar:
    cb=plt.colorbar(c)
    if color_by=='r':
        cbarlabel='Distance (km)'
    elif color_by=='mw':
        cbarlabel='M'
    cb.set_label(cbarlabel)
    
    #Set axis limits and labels:
    plt.xlim(axlims[0])
    plt.ylim(axlims[1])
    
    plt.xlabel(termname)
    plt.ylabel(indname)
    plt.title(ptitle)
    
    #Save figure:
    pngname=fig_dir+index+'_'+term+'_'+color_by+'.png'
    pdfname=pdf_dir+index+'_'+term+'_'+color_by+'.pdf'
    plt.savefig(pngname)
    plt.savefig(pdfname)
    
    return f1
    
 
##########
def plot_terms_colored_condition(home,run_name,condition,robj,robj_x,robj_y,term,index,axlims,color_by,color_by_ind,cvals,mymap):
    '''
    Plot path terms vs. some index, with path terms on the x-axis.
    Input:
        home:           String with path to the home directory for project
        run_name:       String with database/inversion combo for residuals
        condition:      String with condition (for plot and filename)
        robj:           Residuals object
        robj_x:         X vector from residuals object
        robj_y:         Y vector from residuals object
        term:           String with term to plot: 'site_terms','W_residual','E_residual'
        index:          STring with index to plot: 'ind'_raytype_modeltype_indextype
                            raytype=p/s; modeltype=vp/vs/vpvs/qp/qs/qpqs;
                            indextype=pathint,normpathint,gradpathint
                            i.e., 'ind_p_vs_normpathint'
        axlims:         Axis limits: [[xmin,xmax],[ymin,ymax]]
        color_by:       String with variable to color by: 'r', or 'mw'
        color_by_ind:   Indices to use for condition, for coloring
        cvals:          Colorbar limits: [cmin,cmax]
        mymap:          String with colormap to use (i.e., 'jet')
    Output:
        Saves plots to home/run_name/figs/pathterm_index
    '''
    from os import path
    import matplotlib.pyplot as plt
    from scipy.stats.stats import pearsonr
    from numpy import array,arange
    import matplotlib.colors as colors
    import matplotlib.cm as cm
    
    #Get run directory, and figure directory:
    run_dir=path.expanduser(home+run_name+'/')
    fig_dir=run_dir+'figs/'
    pdf_dir=fig_dir+'pdfs/'
    
    #Get a string for the plot name and title:
    figbasename=term+index+'_'+color_by+condition
    
    #Index type...
    indtype=index.split('_')[-1]
    if indtype=='pathint':
        indname='Path Integral'
    elif indtype=='normpathint':
        indname='Normalized Path Integral'
    elif indtype=='gradpathint':
        indname='Path Integral of the gradient'
        
    #Term type:
    if term=='site_terms':
        termname='Site Term (ln residual)'
        termtitle='Site Term'
    elif term=='W_residual':
        termname='Within Event residual (ln residual)'
        termtitle='Within-Event Residual'
    elif term=='E_residual':
        termname='Event residual (ln residual)'
        termtitle='Event Residual'
    elif term=='path_terms':
        termname='Path term (ln residuals)'
        termtitle='Path Term'
        
    #What to plot?
    x=robj_x
    y=robj_y
        
    #get correlation coefficient:
    pcoeff,tails=pearsonr(x,y)
    pcoeff=round(pcoeff,2)
    
    #Title:
    ptitle='Plot of '+termtitle+' vs. '+indname+'for '+condition+'\n Pearson coefficient: '+str(pcoeff)
    
    #Get colormap
    #Make colormap:
    colormap_mwdist=plt.get_cmap(mymap)
    #Make a normalized colorscale
    cmin=cvals[0]
    cmax=cvals[1]
    cNorm=colors.Normalize(vmin=cmin, vmax=cmax)
    #Apply normalization to colormap:
    scalarMap=cm.ScalarMappable(norm=cNorm, cmap=colormap_mwdist)
    
    #Make a fake contour plot for the colorbar:
    Z=[[0,0],[0,0]]
    levels=arange(cmin,cmax,0.01)
    c=plt.contourf(Z, levels, cmap=colormap_mwdist)   
    
    #Assign values to colormap
    colorVal = scalarMap.to_rgba(getattr(robj,color_by)[color_by_ind])
    
    #Plot:
    f1=plt.figure()
    plt.scatter(x,y,facecolors='none',edgecolors=colorVal,lw=0.9)
    
    #Add colorbar:
    cb=plt.colorbar(c)
    if color_by=='r':
        cbarlabel='Distance (km)'
    elif color_by=='mw':
        cbarlabel='M'
    cb.set_label(cbarlabel)
    
    #Set axis limits and labels:
    plt.xlim(axlims[0])
    plt.ylim(axlims[1])
    
    plt.xlabel(termname)
    plt.ylabel(indname)
    plt.title(ptitle)
    
    #Save figure:
    pngname=fig_dir+figbasename+'.png'
    pdfname=pdf_dir+figbasename+'.pdf'
    plt.savefig(pngname)
    plt.savefig(pdfname)
    
    return f1
 
 
 
    
       
###########################3             
def grid_path_term(home,run_name,bindims,raytype,stat_type):
    '''
    Average the path terms for every cell on a 3d grid.  Uses run_name_robj_raydat.pckl
    Input:
        home:               String with home directory (i.e., anza/models/residuals)
        run_name:           String with run name, database/inversion combo
        bindims:            Bin dimensions (nx, ny, nz)
        raytype:            Type of ray: 0=Vp, 1=Vs
        stat_type:          String, type of statistic for binning:
                                'mean', 'median', 'count',or sum'
    Output:
        grid_object:            Gridded object, also saves to /home/run_name/run_name_pterm_grid.pckl
    '''
    
    from scipy.stats import binned_statistic_dd 
    import cPickle as pickle
    from numpy import ones,r_,c_
    import cdefs as cdf
    from os import path
    
    run_dir=path.expanduser(home+run_name+'/')
    rpath=run_dir+run_name+'_robj_raydat.pckl'
    gobjpath=run_dir+run_name+'_pterm_grid.pckl'
    
    #Open object:
    rfile=open(rpath,'r')
    robj=pickle.load(rfile)
    rfile.close()
    
    ##Setup input##
    #First make the path term to match the lat and lon format - a list of arrays, with the same path term value:
    path_list=[]
    for ray_i in range(len(robj.path_terms)):
        ray_length_i=len(robj.vs_lon[ray_i])
        path_terms_i=robj.path_terms[ray_i]*ones(ray_length_i)
        
        #Append to the list of path terms:
        path_list.append(path_terms_i)
        
    
    #Now loop over the number of rays, and for each ray, append it to the larger array:
    #Classify if it's p or s:
    #If it's P:
    if raytype==0:
        for ray_i in range(len(path_list)):
            #If it's the first ray, initiate the arrays:
            if ray_i==0:
                lon=robj.vp_lon[ray_i]
                lat=robj.vp_lat[ray_i]
                dep=robj.vp_depth[ray_i]
                path=path_list[ray_i]
            else:
                lon=r_[lon,robj.vp_lon[ray_i]]
                lat=r_[lat,robj.vp_lat[ray_i]]
                dep=r_[dep,robj.vp_depth[ray_i]]
                path=r_[path,path_list[ray_i]]
                
    elif raytype==1:
        for ray_i in range(len(path_list)):
            #If it's the first ray, initiate the arrays:
            if ray_i==0:
                lon=robj.vs_lon[ray_i]
                lat=robj.vs_lat[ray_i]
                dep=robj.vs_depth[ray_i]
                path=path_list[ray_i]
            else:
                lon=r_[lon,robj.vs_lon[ray_i]]
                lat=r_[lat,robj.vs_lat[ray_i]]
                dep=r_[dep,robj.vs_depth[ray_i]]
                path=r_[path,path_list[ray_i]]
                
    #Concatenate width wise to put in:
    sample=c_[lon,lat,dep]
        
    ##Now bin...
    statistic,bin_edges,binnumber=binned_statistic_dd(sample,path,statistic=stat_type,bins=bindims)
    
    #Save as an object:
    grid_object=cdf.pterm_3dgrid(statistic,bin_edges,binnumber)
    #Save:
    gfile=open(gobjpath,'w')
    pickle.dump(grid_object,gfile)
    gfile.close()
    
    #And return:
    return grid_object
    
    
#####################################################################
# Initiate directory for comparing residuals...
def compare_init(home,run_name):
    '''
    Initiate directory for comparing residuals
    Input:
        home:       String with name to project residuals directory (i.e., '/home/vsahakian/anza/models/residuals')
        run_name:   String with name for comparison directory ('anza_5sta_compare)
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
        makedirs(run_dir+'figs/')
        makedirs(run_dir+'figs/pdfs/')
        
        #If it's ok to go ahead and clobber, set the run variable to 1:
        runall=1
        
    if clob is 'n' or clob is 'N':
        #Do not run!!!!
        #Set the run variable to 0:
        runall=0
        
        
    return runall
    
##############################################################################################
def get_comparison_indices(larger_db,smaller_db):
    '''
    VJS 1/17
    Get the indices in the larger database that correspond to the smaller database. 
    Input:
        larger_db:             Database or residual object to compare with more recordings
        smaller_db:            Databse or residual object to compare with fewer recordings
    Output:
        comparison_indices:    Array with indices of the larger database that match the property in teh smaller
    '''
    
    import numpy as np
    
    # Initiate comparison index array as a list:
    comparison_indices=[]
    
    # Find where in the larger database do the event/site combos match the smaller:
    for recording_i in range(len(smaller_db.evnum)):
        recording_index=np.where((larger_db.evnum==smaller_db.evnum[recording_i]) & (larger_db.sta==smaller_db.sta[recording_i]))[0]
        # Append to the comparision index list:
        comparison_indices.append(recording_index)
        
    # Turn comparison_indices into an array
    comparison_indices=np.array(comparison_indices)
    
    # Return:
    return comparison_indices
    
    
    
##############################################################################
def num_of_randeffects(robj,numev=None,numsta=None):
    '''
    Get the number of events recorded on each station, or number of stations recording each event
    VJS 1/17
    Input:
        robj:               Database or residuals object
        numev:              If True, gets the number of events recorded at each station
        numsta:             If True, gets the number of statiosn recording each event
    Return:
        num_random_effects: Array with requested number of random effects, specified by numev or numsta     
    '''
    import numpy as np
    
    if numsta==True:
            
        #### How many stations record each unique event? ###
        # Get an array with the number of stations recording each event:
        event_numstas = []
        unique_events,unevind=np.unique(robj.evnum,return_index=True)
        
        for eventi in range(len(unique_events)):
            evwhere = np.where(robj.evnum==unique_events[eventi])[0]
            numstas = len(evwhere)
            # append the number of statoin srecording this event to the list:
            event_numstas.append(numstas)
            
        # At the end, now turn into an array:
        event_numstas = np.array(event_numstas)
        
        return event_numstas
    
    elif numev==True:
        ################
        #### How many events are recorded at each station? ###
        station_numevs = []
        unique_stas,unstind=np.unique(robj.sta,return_index=True)
        
        for stai in range(len(unique_stas)):
            stwhere = np.where(robj.sta==unique_stas[stai])[0]
            numevs = len(stwhere)
            # append the number of events recorded on this station to the list:
            station_numevs.append(numevs)
            
        # Turn into an array:
        station_numevs = np.array(station_numevs)
        
        return station_numevs


##############################################################################
def compare_mixed_traditional(home,run_name,tradpath,mixedpath,mymap,evaxlim,staaxlim,pathaxlim,symbol_size,cbins):
    #VJS 1/17
    
    '''
    Makes plots to compare residuals for the same database using both mixed and traditional methods 
    Input:
        home:           String with path to the home directory for project
        run_name:       String with name for model comparison directory
        axlims:         Axis limits: [[xmin,xmax],[ymin,ymax]]
        color_by:       String with variable to color by: 'r', or 'mw'
        cvals:          Colorbar limits: [cmin,cmax]
        mymap:          String with colormap to use (i.e., 'jet')
        symbol_size:    Array/list with symbol sizes for plots: [event,path,site]
        cbins:          Array/list with bin size for colorscale: [event, site]
    Output:
        Saves plots to home/run_name/figs/pathterm_index
        
    '''
    import matplotlib.pyplot as plt
    import cPickle as pickle
    import numpy as np
    import matplotlib.colors as colors
    import matplotlib.cm as cm
    from os import path    

    
    #Get run directory, and figure directory:
    run_dir=path.expanduser(home+run_name+'/')
    fig_dir=run_dir+'figs/'
    pdf_dir=fig_dir+'pdfs/'
    
    
    # Import models:
    tradfile=open(tradpath,'r')
    trad = pickle.load(tradfile)
    tradfile.close()
    
    mixedfile=open(mixedpath,'r')
    mixed = pickle.load(mixedfile)
    mixedfile.close()
    
    # Get info to plot:
    evnums,evind = np.unique(trad.evnum,return_index=True)
    tradevent = np.array(trad.E_residual)[evind]
    mixedevent = np.array(mixed.E_residual)[evind]
    
    tradpath = trad.path_terms
    mixedpath = mixed.path_terms
    
    sta,stind = np.unique(trad.sta,return_index=True)
    stnum=trad.stnum[stind]
    tradsite = np.array(trad.site_terms)[stind]
    mixedsite = np.array(mixed.site_terms)[stind]
    
    
    #################################
    ######  Auxilliary info    ######
    #################$###############
    
    #### How many stations record each unique event? ###
    # Get an array with the number of stations recording each event:
    event_numstas = []
    unique_events,unevind=np.unique(trad.evnum,return_index=True)
    
    for eventi in range(len(unique_events)):
        evwhere = np.where(trad.evnum==unique_events[eventi])[0]
        numstas = len(evwhere)
        # append the number of statoin srecording this event to the list:
        event_numstas.append(numstas)
        
    # At the end, now turn into an array:
    event_numstas = np.array(event_numstas)
    
    ################
    #### How many events are recorded at each station? ###
    station_numevs = []
    unique_stas,unstind=np.unique(trad.sta,return_index=True)
    
    for stai in range(len(unique_stas)):
        stwhere = np.where(trad.sta==unique_stas[stai])[0]
        numevs = len(stwhere)
        # append the number of events recorded on this station to the list:
        station_numevs.append(numevs)
        
    # Turn into an array:
    station_numevs = np.array(station_numevs)
    
    ### Standard Deviations for plot info:
    trad_path_std = np.std(tradpath)
    mixed_path_std = np.std(mixedpath)
    
    trad_event_std = np.std(tradevent)
    mixed_event_std = np.std(mixedevent)
    
    trad_site_std = np.std(tradsite)
    mixed_site_std = np.std(mixedsite)
    
    
    ##########################################################################################
    ##########################
    ######  Plotting    ######
    #################$########
    
    print 'Plotting event terms...'
    
    ## First event terms ##
    # Make a straight line for relationship:
    xe=np.linspace(evaxlim[0][0],evaxlim[0][1],2)
    ye=np.linspace(evaxlim[1][0],evaxlim[1][1],2)
    
    # Plotting colored by number of stations recording each event:
    cmin = min(event_numstas)
    cmax = max(event_numstas)-4
    
    ##Plot:
    #Get colormap
    #Make colormap:
    colormap_numstas=plt.get_cmap(mymap)
    #Make a normalized colorscale
    cNorm=colors.Normalize(vmin=cmin, vmax=cmax)
    #Apply normalization to colormap:
    scalarMap=cm.ScalarMappable(norm=cNorm, cmap=colormap_numstas)
    
    #Make a fake contour plot for the colorbar:
    Z=[[0,0],[0,0]]
    levels=np.arange(cmin,cmax,cbins[0])
    c=plt.contourf(Z, levels, cmap=colormap_numstas)
    
    #Assign values to colormap
    colorVal = scalarMap.to_rgba(event_numstas)
    
    eventfig = plt.figure()
    #plt.scatter(tradevent,mixedevent,marker='o',s=7,color='#333333')
    plt.scatter(tradevent,mixedevent,marker='o',s=symbol_size[0],color=colorVal)
    plt.plot(xe,ye,color='gray',linewidth=1.5)
    
    # Colrobar:
    cb = plt.colorbar(c)
    cb.set_label('Number of stations per event')
    
    # Axis limits:
    plt.xlim(evaxlim[0][0],evaxlim[0][1])
    plt.ylim(evaxlim[1][0],evaxlim[1][1])
    
    plt.xlabel('Single-mean event term (ln residual) \n Std. Dev: %.2f' % trad_event_std)
    plt.ylabel('Mixed effects inversion event term (ln residual) \n Std. Dev: %.2f' % mixed_event_std)
    plt.title('Plot of Single-mean vs. Mixed effects event terms')
    
    #Show the figures
    eventfig.show()
    
    evpngfile = fig_dir+run_name+'_event_comp.png'
    evpdffile = pdf_dir+run_name+'_event_comp.pdf'
    
    print 'Saving event figures'
    eventfig.savefig(evpngfile)
    eventfig.savefig(evpdffile)
    
    
    
    ###############################################
    ## Then make path terms ##
    xp=np.linspace(pathaxlim[0][0],pathaxlim[0][1],2)
    yp=np.linspace(pathaxlim[1][0],pathaxlim[1][1],2)
    
    print 'Plotting path terms'
    pathfig = plt.figure()
    plt.scatter(tradpath,mixedpath,marker='o',s=symbol_size[1],color='#333333')
    plt.plot(xp,yp,color='#9C9C9C',linewidth=1.5)
    
    # Axis limits:
    plt.xlim(pathaxlim[0][0],pathaxlim[0][1])
    plt.ylim(pathaxlim[1][0],pathaxlim[1][1])
    
    plt.xlabel('Single-mean inversion path term (ln residual) \n Std Dev: %.2f' % trad_path_std)
    plt.ylabel('Mixed effects inversion path term (ln residual) \n Std Dev: %.2f' % mixed_path_std)
    plt.title('Plot of Single-mean vs. Mixed effects path terms')
    
    #Show the figure
    pathfig.show()
    
    #Save the figure
    ppngfile = fig_dir+run_name+'path_comp.png'
    ppdffile = pdf_dir+run_name+'path_comp.pdf'
    
    print 'Saving path figures'
    pathfig.savefig(ppngfile)
    pathfig.savefig(ppdffile)
    
    
    ###############################################
    
    ## Then site terms: ##
    xs=np.linspace(staaxlim[0][0],staaxlim[0][1],2)
    ys=np.linspace(staaxlim[1][0],staaxlim[1][1],2)
    
    # Plotting colored by number of events recorded at each station:
    cmin = min(station_numevs)
    cmax = max(station_numevs)
    
    print 'Plotting the site figure'
    
    print 'Making colormap'
    ##Plot:
    #Get colormap
    #Make colormap:
    colormap_numevs=plt.get_cmap(mymap)
    #Make a normalized colorscale
    cNorm=colors.Normalize(vmin=cmin, vmax=cmax)
    #Apply normalization to colormap:
    scalarMap=cm.ScalarMappable(norm=cNorm, cmap=colormap_numevs)
    
    print 'Making fake figure for colormap'
    #Make a fake contour plot for the colorbar:
    Z=[[0,0],[0,0]]
    levels=np.arange(cmin,cmax,cbins[1])
    c=plt.contourf(Z, levels, cmap=colormap_numevs)
    
    print 'Assigning values to colormap'
    #Assign values to colormap
    colorVal = scalarMap.to_rgba(station_numevs)
    
    print 'starting figure, plotting...'
    stafig = plt.figure()
    #plt.scatter(tradevent,mixedevent,marker='o',s=7,color='#333333')
    plt.scatter(tradsite,mixedsite,marker='o',s=symbol_size[2],color=colorVal)
    plt.plot(xs,ys,color='gray',linewidth=1.5)
    
    print 'Making colorbar'
    # Colrobar:
    cb = plt.colorbar(c)
    cb.set_label('Number of events per station')
    
    print 'Axis limits...'
    
    # Axis limits:
    plt.xlim(staaxlim[0][0],staaxlim[0][1])
    plt.ylim(staaxlim[1][0],staaxlim[1][1])
    
    print 'Labels'
    plt.xlabel('Single-mean inversion site term (ln residual) \n Std Dev: %.2f' % trad_site_std)
    plt.ylabel('Mixed effects inversion site term (ln residual) \n Std Dev: %.2f' % mixed_site_std)
    plt.title('Plot of Single-mean vs. Mixed effects site terms')
    
    stafig.show()
    
    stpngfile = fig_dir+'station_comp.png'
    stpdffile = pdf_dir+run_name+'station_comp.pdf'
    
    print 'Saving station figures'
    stafig.savefig(stpngfile)
    stafig.savefig(stpdffile)
    
    
##################################################################################
def compare_models_2dhist(residuals,num_randeffect,nbins,pltaxis,mymap,axlims,clims,xlabel_toggle,ylabel_toggle):
        '''
        
        Input:
            residuals:      Array with the residuals (unique event terms, site terms, etc.)
            num_randeffect: Array with the number of stations recording each event, or number of events recorded on each station
            nbins:          [number of bins for residuals, number of bins for random effect]
            pltaxis:        Axes to plot on. Can use this to plot iteratively
            mymap:          Colormap RGB values to plot arr([r,g,b])
            axlims:         Array with axis limits: [[xmin, xmax],[ymin,ymax]]
            clims:          Limits for colorbar [cmin, cmax]
            xlabel_toggle:  String with instructions to add x label: 'event'/'site'
            ylabel_toggle:  String with instructions to add y label: 'sta_per_ev'/'ev_per_sta'
        Returns:
            
        '''
        
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt
        import numpy as np
        
        
        
        pltaxis.hist2d(residuals, num_randeffect, nbins, range=axlims, cmin=clims[0],cmax=clims[1],cmap=mymap)
        #pltaxis.colorbar()
        
        
        #  If it's plotting event terms:
        if xlabel_toggle=='event':
            pltaxis.set_xlabel('Event Term - ln Residuals')
        elif xlabel_toggle=='site':
            pltaxis.set_xlabel('Site Term - ln Residuals')
        
        if ylabel_toggle=='sta_per_ev':
            pltaxis.set_ylabel('Number of Stations per Event')
        elif ylabel_toggle=='ev_per_sta':
            pltaxis.set_ylabel('Number of Events per Station')



###########
def z_test(sample1,sample2):
    #VJS 1/2017
    '''
    Get a z-value showing if the two samples/distributions are the same distribution or not.  
    If the Z-statistic is less than 2, the two samples are the same. 
    If the Z-statistic is between 2.0 and 2.5, the two samples are marginally different 
    If the Z-statistic is between 2.5 and 3.0, the two samples are significantly different 
    If the Z-statistic is more then 3.0, the two samples are highly signficantly different
    
    Input:
        sample1:        Array with values of the first sample from distribution1
        sample2:        Array with values of the second sample from distribution2
    Returns:
        z_value:        The z-statistic to show how similar distributions are (read above)
    '''
    
    from numpy import mean, std, sqrt
    mean1 = mean(sample1)
    mean2 = mean(sample2)
    
    std_norm1 = std(sample1)/sqrt(len(sample1))
    std_norm2 = std(sample2)/sqrt(len(sample2))
    
    z_value = (mean1 - mean2)/(sqrt(std_norm1**2 + std_norm2**2))
    
    return z_value
    
    

##############################################################################################
def plot_method_method(home,run_name,method_type,largerpath,smallerpath,largername,smallername,comp_indices,mymap,evaxlim,staaxlim,pathaxlim,symbol_size,cbins,clims):
    '''
    VJS 1/2017
    Plot the single-mean vs. single-mean, or mixed vs. mixed terms (event,site,path) for two databases
    Input:
        home:           String with path to project home
        run_name:       String with the comparison directory name
        method_type:    String with the type of method plotting: 'single-mean'/'mixed'
        largerpath:     String with the path to the database/res object with more recordings
        smallerpath:    String with the path to the database/res object with fewer recordings
        largername:     String with the name of the larger db for plotting
        smallername:    String with the name of the smaller db for plotting
        comp_indices:   Indices for matching the db's (see the function ra.get_comparison_indices)
        mymap:          String with the colormap
        evaxlim:        Axis limits for the event plot [[xmin,xmax],[ymin,ymax]]
        staaxlim:       Axis limits for the station plot[[xmin,xmax],[ymin,ymax]]
        pathaxlim:      Axis limits for the path plot [[xmin,xmax],[ymin,ymax]]
        symbolsize:     Array/list with symbol sizes for plots: [event,path,site]
        cbins:          Array/list with bin size for colorscale: [event, site]
        clims:          Limits for colorbar [[cmin_event,cmax_event],[cmin_sta,cmax_sta]]
    Returns:
        
    '''
        
    import res_analysis as ra
    import numpy as np
    from os import path
    import matplotlib.pyplot as plt
    import cPickle as pickle
    import matplotlib.colors as colors
    import matplotlib.cm as cm

    
    # 3.  event
    # 4.  path...
    # 5.  station...
     
    #Get run directory, and figure directory:
    run_dir=path.expanduser(home+run_name+'/')
    fig_dir=run_dir+'figs/'
    pdf_dir=fig_dir+'pdfs/'
    
    ## Open objects:
    largerfile=open(largerpath,'r')
    lobj=pickle.load(largerfile)
    largerfile.close()
    
    smallerfile=open(smallerpath,'r')
    sobj=pickle.load(smallerfile)
    smallerfile.close()
    
    # Color by the smaller database. Get number of events perstation, and number of stations per event for color...
    numev_persta_smaller=ra.num_of_randeffects(sobj,numev=True)
    numsta_perev_smaller=ra.num_of_randeffects(sobj,numsta=True)
    
    #### Get the stuff to plot ####
    # Event stuff:
    ev_unique,ev_uniq_ind = np.unique(lobj.evnum[comp_indices],return_index=True)
    larger_event = np.array(lobj.E_residual)[comp_indices][ev_uniq_ind]
    smaller_event = np.array(sobj.E_residual)[ev_uniq_ind]
    
    larger_ev_std = np.std(larger_event)
    smaller_ev_std = np.std(smaller_event)
    
    # Path stuff:
    larger_path = np.array(lobj.path_terms)[comp_indices]
    smaller_path = np.array(sobj.path_terms)
    
    larger_pa_std = np.std(larger_path)
    smaller_pa_std = np.std(smaller_path)
    
    # Site stuff:
    st_unique,st_uniq_ind = np.unique(lobj.sta[comp_indices],return_index=True)
    larger_site = np.array(lobj.site_terms)[comp_indices][st_uniq_ind]
    smaller_site = np.array(sobj.site_terms)[st_uniq_ind]
    
    larger_st_std = np.std(larger_site)
    smaller_st_std = np.std(smaller_site)
    
    
    ####################    
    ######  Event ######
    ####################
    
    print 'Plotting event terms...'
    
    ## First event terms ##
    # Make a straight line for relationship:
    xe=np.linspace(evaxlim[0][0],evaxlim[0][1],2)
    ye=np.linspace(evaxlim[1][0],evaxlim[1][1],2)
    
    # Plotting colored by number of stations recording each event:
    cmin = clims[0][0]
    cmax = clims[0][1]
    
    ##Plot:
    #Get colormap
    #Make colormap:
    colormap_numstas=plt.get_cmap(mymap)
    #Make a normalized colorscale
    cNorm=colors.Normalize(vmin=cmin, vmax=cmax)
    #Apply normalization to colormap:
    scalarMap=cm.ScalarMappable(norm=cNorm, cmap=colormap_numstas)
    
    #Make a fake contour plot for the colorbar:
    Z=[[0,0],[0,0]]
    levels=np.arange(cmin,cmax,cbins[0])
    c=plt.contourf(Z, levels, cmap=colormap_numstas)
    
    #Assign values to colormap
    colorVal = scalarMap.to_rgba(numsta_perev_smaller)
    print 'max of numsta_perev is %f' % max(numsta_perev_smaller)
    print 'min of numsta_perev is %f' % min(numsta_perev_smaller)
    
    eventfig = plt.figure()
    #plt.scatter(tradevent,mixedevent,marker='o',s=7,color='#333333')
    plt.scatter(smaller_event,larger_event,marker='o',s=symbol_size[0],color=colorVal)
    plt.plot(xe,ye,color='gray',linewidth=1.5)
    
    # Colrobar:
    cb = plt.colorbar(c)
    cb.set_label('Number of stations per event')
    
    # Axis limits:
    plt.xlim(evaxlim[0][0],evaxlim[0][1])
    plt.ylim(evaxlim[1][0],evaxlim[1][1])
    
    if method_type=='single-mean':
        plt.xlabel('Single-mean event term, %s (ln residual) \n Std. Dev: %.2f' % (smallername,smaller_ev_std))
        plt.ylabel('Single-mean event term, %s (ln residual) \n Std. Dev: %.2f' % (largername,larger_ev_std))
        plt.title('Plot of Single-mean vs. Single-mean event terms')
    
        evpngfile = fig_dir+run_name+'_'+smallername+'_'+largername+'_event_single_comp.png'
        evpdffile = pdf_dir+run_name+'_'+smallername+'_'+largername+'_event_single_comp.pdf'
        print 'Saving figure to %s' % evpngfile
        
    elif method_type=='mixed':
        plt.xlabel('Mixed event term, %s (ln residual) \n Std. Dev: %.2f' % (smallername,smaller_ev_std))
        plt.ylabel('Mixed event term, %s (ln residual) \n Std. Dev: %.2f' % (largername,larger_ev_std))
        plt.title('Plot of Mixed vs. Mixed event terms')
    
        evpngfile = fig_dir+run_name+'_'+smallername+'_'+largername+'_event_mixed_comp.png'
        evpdffile = pdf_dir+run_name+'_'+smallername+'_'+largername+'_event_mixed_comp.pdf'
        print 'Saving figure to %s' % evpngfile
    
    print 'Saving event figures'
    eventfig.savefig(evpngfile)
    eventfig.savefig(evpdffile)
    
    #Show the figures
    eventfig.show()
    
    
    ####################
    ###### Station #####
    ####################
    
    print 'Plotting station terms...'
    
    ## First event terms ##
    # Make a straight line for relationship:
    xs=np.linspace(staaxlim[0][0],staaxlim[0][1],2)
    ys=np.linspace(staaxlim[1][0],staaxlim[1][1],2)
    
    # Plotting colored by number of stations recording each event:
    cmin = clims[1][0]
    cmax = clims[1][1]
    
    ##Plot:
    #Get colormap
    #Make colormap:
    colormap_numevs=plt.get_cmap(mymap)
    #Make a normalized colorscale
    cNorm=colors.Normalize(vmin=cmin, vmax=cmax)
    #Apply normalization to colormap:
    scalarMap=cm.ScalarMappable(norm=cNorm, cmap=colormap_numevs)
    
    #Make a fake contour plot for the colorbar:
    Z=[[0,0],[0,0]]
    levels=np.arange(cmin,cmax,cbins[1])
    c=plt.contourf(Z, levels, cmap=colormap_numevs)
    
    #Assign values to colormap
    colorVal = scalarMap.to_rgba(numev_persta_smaller)
    print 'max of numev_persta is %f' % max(numev_persta_smaller)
    print 'min of numev_persta is %f' % min(numev_persta_smaller)
    
    stafig = plt.figure()
    #plt.scatter(tradevent,mixedevent,marker='o',s=7,color='#333333')
    plt.scatter(smaller_site,larger_site,marker='o',s=symbol_size[0],color=colorVal)
    plt.plot(xs,ys,color='gray',linewidth=1.5)
    
    # Colrobar:
    cb = plt.colorbar(c)
    cb.set_label('Number of events per station')
    
    # Axis limits:
    plt.xlim(staaxlim[0][0],staaxlim[0][1])
    plt.ylim(staaxlim[1][0],staaxlim[1][1])
    
    if method_type=='single-mean':
        plt.xlabel('Single-mean site term, %s (ln residual) \n Std. Dev: %.2f' % (smallername,smaller_st_std))
        plt.ylabel('Single-mean site term, %s (ln residual) \n Std. Dev: %.2f' % (largername,larger_st_std))
        plt.title('Plot of Single-mean vs. Single-mean site terms')
    
        stpngfile = fig_dir+run_name+'_'+smallername+'_'+largername+'_site_single_comp.png'
        stpdffile = pdf_dir+run_name+'_'+smallername+'_'+largername+'_site_single_comp.pdf'
        print 'Saving figure to %s' % stpngfile
        
    elif method_type=='mixed':
        plt.xlabel('Mixed site term, %s (ln residual) \n Std. Dev: %.2f' % (smallername,smaller_st_std))
        plt.ylabel('Mixed site term, %s (ln residual) \n Std. Dev: %.2f' % (largername,larger_st_std))
        plt.title('Plot of Mixed vs. Mixed site terms')
    
        stpngfile = fig_dir+run_name+'_'+smallername+'_'+largername+'_site_mixed_comp.png'
        stpdffile = pdf_dir+run_name+'_'+smallername+'_'+largername+'_site_mixed_comp.pdf'
        print 'Saving figure to %s' % stpngfile
        
    print 'Saving site figures'
    stafig.savefig(stpngfile)
    stafig.savefig(stpdffile)
    
    #Show the figures
    stafig.show()
    
    
    
    ####################
    ####### Path #######
    ####################
    
    print 'Plotting path terms...'
    
    ## First event terms ##
    # Make a straight line for relationship:
    xp=np.linspace(pathaxlim[0][0],pathaxlim[0][1],2)
    yp=np.linspace(pathaxlim[1][0],pathaxlim[1][1],2)
    
    # Path
    pafig = plt.figure()
    #plt.scatter(tradevent,mixedevent,marker='o',s=7,color='#333333')
    plt.scatter(smaller_path,larger_path,marker='o',s=symbol_size[0],color='#333333')
    plt.plot(xp,yp,color='gray',linewidth=1.5)
    
    
    # Axis limits:
    plt.xlim(pathaxlim[0][0],pathaxlim[0][1])
    plt.ylim(pathaxlim[1][0],pathaxlim[1][1])
    
    if method_type=='single-mean':
        plt.xlabel('Single-mean path term, %s (ln residual) \n Std. Dev: %.2f' % (smallername,smaller_pa_std))
        plt.ylabel('Single-mean path term, %s (ln residual) \n Std. Dev: %.2f' % (largername,larger_pa_std))
        plt.title('Plot of Single-mean vs. Single-mean path terms')
    
        papngfile = fig_dir+run_name+'_'+smallername+'_'+largername+'_path_single_comp.png'
        papdffile = pdf_dir+run_name+'_'+smallername+'_'+largername+'_path_single_comp.pdf'
        print 'Saving figure to %s' % papngfile
        
    elif method_type=='mixed':
        plt.xlabel('Mixed path term, %s (ln residual) \n Std. Dev: %.2f' % (smallername,smaller_pa_std))
        plt.ylabel('Mixed path term, %s (ln residual) \n Std. Dev: %.2f' % (largername,larger_pa_std))
        plt.title('Plot of Mixed vs. Mixed path terms')
        
        papngfile = fig_dir+run_name+'_'+smallername+'_'+largername+'_path_mixed_comp.png'
        papdffile = pdf_dir+run_name+'_'+smallername+'_'+largername+'_path_mixed_comp.pdf'
        print 'Saving figure to %s' % papngfile
    
    print 'Saving path figures'
    pafig.savefig(papngfile)
    pafig.savefig(papdffile)
    
    #Show the figures
    pafig.show()


#############################################################################
def plot_binned_metric(residualobj,residualterm,metric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz):
    '''
    VJS 10/2017
    Plots a metric vs. residual for certain bins, returns a figure with number
    of subplots equal to number of bins
    Input: 
        residualobj:            Residuals object, including metrics
        residualterm:           String with residual to plot: Event - 'event', Path - 'path', Site - 'site'
        metric:                 String with metric to plot: Path integral - 'ind_{ray - p or s}_{velocity - vp or vs}_pathint', Normalized path integral - 'normpathint', Integral of gradient - 'gradpathint'
        binedges:               Array with bin edges: ([bin1left,bin1right,bin2right,bin3right,etc.])
        bin_by:                 String with variable to bin by: Rrup - 'Rrup', M - 'M', Az from site - 'Site Azimuth'
        axlims:                 List with axis limits: [[xmin,xmax],[ymin,ymax]]
        color_by:               String with variable to color by: Rrup - 'Rrup', M - 'M', Az from site - 'Site Azimuth'
        colorscheme:            String with colormap to use 
        clims:                  List with colorscale limits: [cmin,cmax,cincrement]
        cbartick:               Tick increment for colorbar (a single number)
        plotdims:               Plot dimensions: [width,height]
        plotrowscols:           List with rows and columns: [numberrows, numbercolumns]
        fontsz:               Number with font size, i.e., 14
    Output:
        binned_figure:          Figure with subplots for bins
    '''
    
    import matplotlib.pyplot as plt
    from numpy import where,arange,around
    from pyproj import Geod
    from scipy.stats.stats import pearsonr
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    
    # Number of bins:
    numberbins = len(binedges)-1
    
    # Get what you're binning by:
    if bin_by == 'Rrup':
        binvalue = residualobj.r
    elif bin_by == 'M':
        binvalue = residualobj.mw
    elif bin_by == 'Site Azimuth':
        # Initiate projection object:
        g = Geod(ellps='WGS84')
        az,backaz,distance = g.inv(residualobj.stlon,residualobj.stlat,residualobj.elon,residualobj.elat)
        binvalue = az
        
    # Get what you're coloring by:
    if color_by == 'Rrup':
        colorvalue = residualobj.r
    elif color_by == 'M':
        colorvalue = residualobj.mw
    elif color_by == 'Site Azimuth':
        # Initiate projection object:
        g = Geod(ellps='WGS84')
        az,backaz,distance = g.inv(residualobj.stlon,residualobj.stlat,residualobj.elon,residualobj.elat)
        colorvalue = az
        
    # Get the metric:
    if metric == 'ind_s_vs_pathint':
        metricvalue = residualobj.ind_s_vs_pathint
        xlabel = 'Path Integral Metric'
    elif metric == 'ind_s_vs_normpathint':
        metricvalue = residualobj.ind_s_vs_normpathint
        xlabel = 'Normalized Path Metric'
    elif metric == 'ind_s_vs_gradpathint':
        metricvalue = residualobj.ind_s_vs_gradpathint
        xlabel = 'Gradient Metric'

    elif metric == 'ind_p_vs_pathint':
        metricvalue = residualobj.ind_p_vs_pathint
        xlabel = 'Path Integral Metric'
    elif metric == 'ind_p_vs_normpathint':
        metricvalue = residualobj.ind_p_vs_normpathint
        xlabel = 'Normalized Path Metric'
    elif metric == 'ind_p_vs_gradpathint':
        metricvalue = residualobj.ind_p_vs_gradpathint
        xlabel = 'Gradient Metric'
        
    # Get the residual:
    if residualterm == 'event':
        residual = residualobj.E_residual
        ylabel = 'Event Residual'
    elif residualterm == 'site':
        residual = residualobj.site_terms
        ylabel = 'Site Residual'
    elif residualterm == 'path':
        residual = residualobj.path_terms
        ylabel = 'Path Residual'
    
    # Now go through bins, and append them to a list:
    residual_list = []
    metric_list = []
    r_value_list = []
    p_value_list = []
    color_value_list = []
    
    # Loop through bins to get indices for each bin:
    for bin_i in range(numberbins):
        i_bin_ind = where((binvalue > binedges[bin_i]) & (binvalue <= binedges[bin_i+1]))
        i_residual = residual[i_bin_ind]
        i_metric = metricvalue[i_bin_ind]
        i_colorval = colorvalue[i_bin_ind]
        
        # Append to their lists:
        residual_list.append(i_residual)
        metric_list.append(i_metric)
        color_value_list.append(i_colorval)
        
        # Get Statistics - pearsons and p value:
        i_bin_r,i_bin_p = pearsonr(i_residual,i_metric)
        
        # Append statistics to lists:
        r_value_list.append(i_bin_r)
        p_value_list.append(i_bin_p)
    
    
    ########################
    
    ## Initiate plot:
    binnedplot, binnedaxes = plt.subplots(nrows=plotrowscols[0],ncols=plotrowscols[1],figsize=(plotdims[0],plotdims[1]))
    
    # Flatten the axes array so that it can easily be looped over per bin:
    binnedaxes = binnedaxes.flatten()
    
    # For each bin, plot:
    for subplot_i in range(numberbins):
        axis_i = binnedaxes[subplot_i]
        
        plottext = str(binedges[subplot_i]) + ' < ' + bin_by + ' <= ' + str(binedges[subplot_i + 1])
        p_val = '%.1e' % p_value_list[subplot_i]
        plottextstats = 'r = ' + str(around(r_value_list[subplot_i],1)) + ', p = ' + str(p_val)
        
        # Make colormap
        colormap=plt.get_cmap(colorscheme)
        # Make a normalized colorscale
        cNorm=colors.Normalize(vmin=clims[0], vmax=clims[1])
        # Apply normalization to colormap:
        scalarMap=cm.ScalarMappable(norm=cNorm, cmap=colormap)
        
        # Make a fake contour plot for the colorbar:
        Z=[[0,0],[0,0]]
        levels=arange(clims[0],clims[1],clims[2])
        c=plt.contourf(Z, levels, cmap=colormap)
        
        # Assign values to colormap
        colorVal = scalarMap.to_rgba(color_value_list[subplot_i])
        
        # Scatter:
        axis_i.scatter(metric_list[subplot_i],residual_list[subplot_i],facecolors='none',edgecolors=colorVal,marker='o',s=5)
        
        # Colorbar: 
        axis_i_cb = plt.colorbar(c,ax=axis_i,shrink=0.70)
        axis_i_cb.set_label(color_by,labelpad=0,rotation=90)
        axis_i_cb.set_ticks(ticks=arange(clims[0],clims[1],cbartick))
        
        # Labels:
        axis_i.set_xlabel(xlabel,labelpad=0,fontsize=fontsz)
        axis_i.set_ylabel(ylabel,labelpad=0,fontsize=fontsz)
        axis_i.set_title(plottext + '\n' + plottextstats,fontsize=fontsz+2)        

        
    # Remove empty axes:
    for axis_j in range(len(binnedaxes)):
        if axis_j >= numberbins:
            binnedaxes[axis_j].axis('off')
            binnedplot.delaxes(binnedaxes[axis_j])
            
            
    # Return figure:
    return binnedplot