##Run Residual analysis

#1 - Make material model object (read data, make, save)
def make_material_object(coordspath,materialpath,matobjectpath,lon_convert):
    '''
    Make a material model object given coordinates and material file; Writes object to file.
    Input:
        coordspath:             String to path of the coordinates file
        materialpath:           String to tpath of the material model file
        matobjectpath:          String to path of output material object location
        lon_convert:            Convert longitude?
                                0=none, 1=-west to +west, 2=+west to -west
    Output:
        materialobject:         Material model object
    '''
    
    import dread as dr
    import cdefs as cdf
    import cPickle as pickle
    
    #Read in data:
    x,y,z,nx,ny,nz,mod=dr.read_material_model(coordspath,materialpath)
    
    #Convert longitude:
    if lon_convert==1:     # negative west to positive west
        #plus x becuase x is negative
        x=360+x
    elif lon_convert==2:   # positive west to negative west
        x=x-360

    #Make object:
    materialobject=cdf.material_model(x,y,z,nx,ny,nz,mod)
    
    #Save it to a file:
    mfile=open(matobjectpath,'w')
    pickle.dump(materialobject,mfile)
    mfile.close()

    #Return object:
    return materialobject


#2 - Interpolate ray values through model
def interp_rays(residualobject,materialobject,interptype,raytypes):
    '''
    Interpolate P and S wave rays through the material model
    Input:
        residualobject:         Object with the residuals, ray positions, etc.
        materialobject:         Object with the materialmodel
        interptype:             String with type of interpolation for Rbf (i.e., 'linear')
        raytypes:               Array with ray types to interpolate: ([0]), or ([0,1]) - 0=Vp, 1=Vs
    Output:
        pvals:                  List of arrays with values of model for P wave rays
        svals:                  List of arrays with values of model for S wave rays
    '''
    
    import res_analysis as ra
    
    #Find which rays to interpolate:
    
    
    #First run for P-waves:
    #INterpolate:
    
    ######   TEST    #########
    #if (len(raytypes)==2) & (0 in raytypes) & (1 in raytypes):
    #    pvals=1
    #    svals=2
    #    
    #    return pvals,svals
    #    
    #elif (len(raytypes)==1) & (0 in raytypes):
    #    pvals=1
    #    
    #    return pvals
    #    
    #elif (len(raytypes)==1) & (1 in raytypes):
    #    svals=2
    #    
    #    return svals
    
    
    
    if (len(raytypes)==2) & (0 in raytypes) & (1 in raytypes):
        rayflag=0
        pvals=ra.interpolate_rays(residualobject,materialobject,interptype,rayflag)
        
        rayflag=1
        svals=ra.interpolate_rays(residualobject,materialobject,interptype,rayflag)
        
        return pvals,svals
        
    elif (len(raytypes)==1) & (0 in raytypes):
        rayflag=0
        pvals=ra.interpolate_rays(residualobject,materialobject,interptype,rayflag)
        
        return pvals
        
    elif (len(raytypes)==1) & (1 in raytypes):
        rayflag=1
        svals=ra.interpolate_rays(residualobject,materialobject,interptype,rayflag)
        
        return svals
    

        

#3 - Compute indices
def compute_index(pvals,svals,materialobject,index_type):
    '''
    Compute indices and save to object
    Input:
        pvals:              List of arrays with values of P rays in material model
        svals:              List of arrays with values of S rays in material model
        materialobject:     Object with material model
        index_type:         Type of index to compute: 0=path integral, 
                                1=normalized path integral, 2=path integral of gradient,
                                3=normalized path integral of gradient
    Output:
        
    '''
    
    

#4 - Save as object