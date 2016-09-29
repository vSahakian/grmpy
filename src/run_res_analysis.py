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
def interp_rays(residualobject,materialobject,interptype):
    '''
    Interpolate P and S wave rays through the material model
    Input:
        residualobject:         Object with the residuals, ray positions, etc.
        materialobject:         Object with the materialmodel
        interptype:             String with type of interpolation for Rbf (i.e., 'linear')
    Output:
        pvals:                  List of arrays with values of model for P wave rays
        svals:                  List of arrays with values of model for S wave rays
    '''
    
    import res_analysis as ra
    
    #First run for P-waves:
    rayflag=0
    
    #INterpolate:
    pvals=ra.interpolate_rays(residualobject,materialobject,interptype,rayflag)
    
    #Now for S-waves:
    rayflag=1
    svals=ra.interpolate_rays(residualobject,materialobject,interptype,rayflag)
    
    #Return
    return pvals, svals
        

#3 - Compute indices
def compute_index(ray_vals,materialobject,index_type):
    '''
    Compute an index


#4 - Save as object