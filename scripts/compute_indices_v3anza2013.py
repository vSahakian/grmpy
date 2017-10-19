# Compute indices for v3anza database
# VJS 10/2017

import cPickle as pickle
import run_res_analysis as runra
import res_analysis as ra

#Change this parameter depending on where you run:
#0=desktop
#1=mac

what_home=1

if what_home==0:
    #Desktop:
    HOME='/media/vsahakian/katmai'
elif what_home==1:
    #Mac:
    HOME='/Users/vsahakian'
    
    
########################
#######  Paths  ########
########################

home=HOME+'/anza/models/residuals/'
run_name='mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30'

## Paths to Fang 2016 model
coordspath=HOME+'/anza/data/vm/Fang2016/coords.txt'
materialmodelpath=HOME+'/anza/data/vm/Fang2016/Vs.dat'

## Path to residuals object:
rpath=HOME+'/anza/models/residuals/mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_raydat.pckl'

## Path to store materialobject:
materialmodelname='FangVs'
# Name of the material model object location
matobjpath=HOME+'/anza/data/pckl/'+materialmodelname+'.pckl'

## Material flag - what kind is it?  0=Vp, 1=Vs, 2=Vp/Vs, 3=Qp, 4=Qs, 5=Qp/Qs
materialflag=1 # Vs
raytype=([0,1]) # Vs

## Residual Object paths:
#interpolated ray values
rpath_interp=rpath.split('.pckl')[0]+'_interp_'+materialmodelname+'.pckl'

#interpolated ray values and indices
rpath_indices=rpath.split('.pckl')[0]+'_interp_indices_'+materialmodelname+'.pckl'



############################################################################

###############################
##### Read in objects #########
###############################

## Read in the residuals object
rfile = open(rpath,'r')
robj = pickle.load(rfile)
rfile.close()


## Read in the material properties object - Fang 2016 Vs
mfile = open(matobjpath,'r')
mobj = pickle.load(mfile)
mfile.close()


###############################
#### Interpolate rays #########
###############################


## Interpolate Vp and Vs rays through this model:
interpolation_type='linear'
p_ray_data,s_ray_data=runra.interp_rays(robj,mobj,interpolation_type,raytype)

## Add to and save in a residuals object:
# Values for p rays:
robj.add_material_values(p_ray_data,materialflag,0)
# Values for s rays:
robj.add_material_values(s_ray_data,materialflag,1)

# Save to file:
rfile_interp=open(rpath_interp,'w')
pickle.dump(robj,rfile_interp)
rfile_interp.close()

###############################
####  Compute Metrics #########
###############################

##Compute indices and save to a new residuals object:
#First compute for the p-waves:
#No normalization:
ind_p_vs_path=ra.compute_pathintegral(p_ray_data,mobj,0)
#normalization:
ind_p_vs_normpath=ra.compute_pathintegral(p_ray_data,mobj,1)
#gradient:
ind_p_vs_gradpath=ra.compute_devpathintegral(p_ray_data,mobj,0)

#Now for s-wves:
ind_s_vs_path=ra.compute_pathintegral(s_ray_data,mobj,0)
#normalization:
ind_s_vs_normpath=ra.compute_pathintegral(s_ray_data,mobj,1)
#gradient:
ind_s_vs_gradpath=ra.compute_devpathintegral(s_ray_data,mobj,0)


##Save the indices:
#indext type =0, ray type is p=0,material type is vs=1) 
robj.add_indices(ind_p_vs_path,0,0,1)
#indextype =1,
robj.add_indices(ind_p_vs_normpath,1,0,1)
#indextype=2:
robj.add_indices(ind_p_vs_gradpath,2,0,1)

#indext type =0, ray type is s=1,material type is vs=1) 
robj.add_indices(ind_s_vs_path,0,1,1)
#indextype =1,
robj.add_indices(ind_s_vs_normpath,1,1,1)
#indextype=2:
robj.add_indices(ind_s_vs_gradpath,2,1,1)


# Save to file:
rfile_indices=open(rpath_indices,'w')
pickle.dump(robj,rfile_indices)
rfile_indices.close()
