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
            
            