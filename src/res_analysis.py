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
    
    tradpath = mixed.path_terms
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
    
    plt.xlabel('Simple inversion event term (ln residual)')
    plt.ylabel('Mixed effects inversion event term (ln residual)')
    plt.title('Plot of Simple vs. Mixed effects event terms')
    
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
    
    plt.xlabel('Simple inversion path term (ln residual)')
    plt.ylabel('Mixed effects inversion path term (ln residual)')
    plt.title('Plot of Simple vs. Mixed effects path terms')
    
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
    plt.xlabel('Simple inversion site term (ln residual)')
    plt.ylabel('Mixed effects inversion site term (ln residual)')
    plt.title('Plot of Simple vs. Mixed effects site terms')
    
    stafig.show()
    
    stpngfile = fig_dir+'station_comp.png'
    stpdffile = pdf_dir+run_name+'station_comp.pdf'
    
    print 'Saving station figures'
    stafig.savefig(stpngfile)
    stafig.savefig(stpdffile)
    
    
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
    