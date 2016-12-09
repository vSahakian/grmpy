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
    
    