######Residual Analysis######
##VJS 6/2016

##Interpolate rays through a material model##
def interpolate_rays(residobj,materialobj,rayflag):
    '''
    Interpolate rays through a material model
    Input:
        residobj:           Object with residuals and ray position information
        materialobj:        Object with material model.  Depth should be positive for input.
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
    interpolator = Rbf(columnx, columny, columnz, data,function='linear')

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
            
            
    #Return info:
    return ray_data
        


      