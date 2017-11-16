###########################################################################
## Compute metrics for set distances away from teh site


import cPickle as pickle
import run_res_analysis as runra
import res_analysis as ra
import numpy as np
import pandas as pd
from pyproj import Geod

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


metricdfpath = home + 'mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30_siteradiusmetrics.pckl'


########################
#######  Params  #######
########################

# BAse model name:
basepath = rpath.split('.pckl')[0]

## Material flag - what kind is it?  0=Vp, 1=Vs, 2=Vp/Vs, 3=Qp, 4=Qs, 5=Qp/Qs
materialflag=1 # Vs

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

# First get azimuths frm site to event:
g = Geod(ellps='WGS84')
azimuth,backazimuth,distance = g.inv(robj.stlon,robj.stlat,robj.elon,robj.elat)

## Make empty lists for metrics:
ind_s_vs_path = np.zeros((len(robj.path_terms),len(site_radius)))
ind_s_vs_normpath = np.zeros((len(robj.path_terms),len(site_radius)))
ind_s_vs_gradpath = np.zeros((len(robj.path_terms),len(site_radius))) 

column_list_path = []
column_list_normpath = []
column_list_gradpath = []


## For every distance metric, get the new rays and a new object:
for radiusi in range(len(site_radius)):
    
    # Add to columns list for dataframe:
    column_list_path.append('path' + np.str(site_radius[radiusi]))
    column_list_normpath.append('normpath' + np.str(site_radius[radiusi]))
    column_list_gradpath.append('gradpath' + np.str(site_radius[radiusi]))
    
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
    ind_s_vs_path[:,radiusi] = i_ind_s_vs_path
    ind_s_vs_normpath[:,radiusi] = i_ind_s_vs_normpath
    ind_s_vs_gradpath[:,radiusi] = i_ind_s_vs_gradpath
    
    ## Save to a file:
    #i_rpath = basepath + '_site' + np.str(site_radius[radiusi]) + 'km.pckl'
    #irfile = open(i_rpath,'w')
    #pickle.dump(i_residobj_radius,irfile)
    #irfile.close()
    
## Now save the metrics themselves to a file 
# First make a dataframe...
#   Start with required terms:
basic_df_array = np.c_[robj.evnum,robj.stnum,robj.mw,robj.r,azimuth,robj.E_residual,robj.site_terms,robj.path_terms,robj.site_stderr]
basic_df = pd.DataFrame(basic_df_array,columns=['evnum','stnum','M','rrup','site2ev_azimuth','E_residual','site_terms','path_terms','site_stderr'])

# Then make array with each metric for all radii:
ind_s_vs_path_df = pd.DataFrame(ind_s_vs_path,columns=column_list_path)
ind_s_vs_normpath_df = pd.DataFrame(ind_s_vs_normpath,columns=column_list_normpath)
ind_s_vs_gradpath_df = pd.DataFrame(ind_s_vs_gradpath,columns=column_list_gradpath)

# Append them all:
site_radius_df = pd.concat([basic_df,ind_s_vs_path_df,ind_s_vs_normpath_df,ind_s_vs_gradpath_df],axis=1)

## Now dump this to an object:
metricfile = open(metricdfpath,'w')
pickle.dump(site_radius_df,metricfile)
metricfile.close()
