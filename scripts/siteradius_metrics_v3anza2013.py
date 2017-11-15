###########################################################################
## Compute metrics for set distances away from teh site
## Compare them with path terms

import cPickle as pickle
import run_res_analysis as runra
import res_analysis as ra
import numpy as np

#Change this parameter depending on where you run:
#0=desktop
#1=mac

what_home=0

if what_home==0:
    #Desktop:
    HOME='/katmai'
elif what_home==1:
    #Mac:
    HOME='/Users/vsahakian'
    
    
########################
#######  Paths  ########
########################

home=HOME+'/anza/models/residuals/'


## Main residual model:
rpath = home + 'mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_raydat.pckl'

## Material model:
mpath = HOME+'/anza/data/pckl/FangVs.pckl'

## Material flag - what kind is it?  0=Vp, 1=Vs, 2=Vp/Vs, 3=Qp, 4=Qs, 5=Qp/Qs
materialflag=1 # Vs


########################
#######  Params  #######
########################

# BAse model name:
basepath = rpath.split('.pckl')[0]

dist_from = 'site'
site_radius = np.array([0.5,1,2,5,10,20])
ray_type='vs'
interpolation_type='linear'
interpraytype = np.array([1])


#################################################################################
#################################################################################

###############################
######  Read in models  #######
###############################

## Read in the residuals object
rfile = open(rpath,'r')
robj = pickle.load(rfile)
rfile.close()
print 'Read in residuals model'


## Read in the material properties object - Fang 2016 Vs
mfile = open(mpath,'r')
mobj = pickle.load(mfile)
mfile.close()
print 'Read in material object model'


###############################
######   Get new rays   #######
###############################

## Make empty lists for metrics:
ind_s_vs_path = []
ind_s_vs_normpath = []
ind_s_vs_gradpathint = [] 

## For every distance metric, get the new rays and a new object:
for radiusi in range(len(site_radius)):

    ## Get new locations and interpolation:
    print 'Getting ray locations for site radius ' + np.str(site_radius[radiusi]) + '...'
    i_radius_distance_lon, i_radius_distance_lat, i_radius_distance_depth, i_residobj_radius = ra.get_rays_inradius(robj,site_radius[radiusi],dist_from,ray_type)
    print 'Have new res object for site radius ' + np.str(site_radius[radiusi]) + '.'
    
    ## INterpolate for metrics...
    print 'About to interpolate rays through velocity model for site radius ' + np.str(site_radius[radiusi])
    i_s_ray_data=runra.interp_rays(i_residobj_radius,mobj,interpolation_type,interpraytype)
    print 'Finished interpolating for site radius '+ np.str(site_radius[radiusi])

    ## Add to and save in a residuals object:
    # Values for s rays:
    print 'Adding to the residual object...'
    i_residobj_radius.add_material_values(i_s_ray_data,materialflag,1)

    ## Compute metrics:
    print 'Computing path integral for site radius ' + np.str(site_radius[radiusi])
    i_ind_s_vs_path=ra.compute_pathintegral(i_s_ray_data,mobj,0)

    #normalization:
    print 'Computing normalized path integral for site radius ' + np.str(site_radius[radiusi])
    i_ind_s_vs_normpath=ra.compute_pathintegral(i_s_ray_data,mobj,1)

    #gradient:
    print 'Computing velocity gradient metric for site radius ' + np.str(site_radius[radiusi])
    i_ind_s_vs_gradpath=ra.compute_devpathintegral(i_s_ray_data,mobj,0)

    ## Save to residuals object:
    #indext type =0, ray type is s=1,material type is vs=1) 
    print 'Saving metrics to residuals object...'
    i_residobj_radius.add_indices(i_ind_s_vs_path,0,1,1)
    #indextype =1,
    i_residobj_radius.add_indices(i_ind_s_vs_normpath,1,1,1)
    #indextype=2:
    i_residobj_radius.add_indices(i_ind_s_vs_gradpath,2,1,1)
    
    # Append to lists:
    ind_s_vs_path.append(i_ind_s_vs_path)
    ind_s_vs_normpath.append(i_ind_s_vs_normpath)
    ind_s_vs_gradpathint.append(i_ind_s_vs_gradpath)
    
    # Save to a file:
    i_rpath = basepath + '_site' + np.str(site_radius[radiusi]) + 'km.pckl'
    irfile = open(i_rpath,'w')
    pickle.dump(i_residobj_radius,irfile)
    irfile.close()
    
#