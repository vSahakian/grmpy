#############  Script to compute Q metrics  ###########
## VJS 2/2018

import dread
import cdefs as cdf
import cPickle as pickle
from numpy import meshgrid,zeros,where,array,genfromtxt
import run_res_analysis as runra
import res_analysis as ra
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr

#######################################
##############  Paths #################
#######################################

#Change this parameter depending on where you run:
#0=desktop
#1=mac

what_home=1

if what_home==0:
    #Desktop:
    HOME='/katmai'
elif what_home==1:
    #Mac:
    HOME='/Users/vsahakian'


## Path to Qp and Qs
hauksson_model_path_Qp = HOME + '/anza/data/vm/Hauksson/Qp_simple.txt'
hauksson_model_path_Qs = HOME + '/anza/data/vm/Hauksson/Qs_simple.txt'

hauksson_grid_path_Qp = HOME + '/anza/data/vm/Hauksson/Qp.pckl'
hauksson_grid_path_Qs = HOME + '/anza/data/vm/Hauksson/Qs.pckl'

#Path to residuals object:
rpath = HOME+'/anza/models/residuals/mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_raydat.pckl'

#Material flag - what kind is it?  0=Vp, 1=Vs, 2=Vp/Vs, 3=Qp, 4=Qs, 5=Qp/Qs
materialflag_Qp = 3
materialflag_Qs = 4

# Material model name
materialmodelname_Qp_Qs = 'HaukQp_Qs'

#Residual Object paths:
#interpolated ray values
rpath_interp_Qp_Qs = rpath.split('.pckl')[0]+'_interp_'+materialmodelname_Qp_Qs+'.pckl'
#rpath_interp_Qp_Qs = '/home/vsahakian/Desktop/' + rpath.split('.pckl')[0].split('residuals/')[1].split('/')[1] +'_interp_'+materialmodelname_Qp_Qs+'.pckl'

#interpolated ray values and indices
rpath_indices_Qp_Qs = rpath.split('.pckl')[0]+'_interp_indices_'+materialmodelname_Qp_Qs+'.pckl'
#tmp_rpath_indices_Qp_Qs = '/home/vsahakian/Desktop/' + rpath.split('.pckl')[0].split('residuals/')[1].split('/')[1] +'_interp_indices_'+materialmodelname_Qp_Qs+'.pckl'


#########################################################
############### First Interpolate rays  #################
#########################################################
#
###Interpolate rays:
##Read in residuals object
#print 'opening residuals object... \n'
#rfile=open(rpath,'r')
#robj=pickle.load(rfile)
#rfile.close()
#print 'residuals object opened \n'
#
#print 'opening Qp... \n'
#materialmodel_Qp = genfromtxt(hauksson_model_path_Qp)
#print 'opened Qp model. \n'
#
#print 'opening Qs... \n'
#materialmodel_Qs = genfromtxt(hauksson_model_path_Qs)
#print 'opened Qs model. \n'
#
#print 'Interpolating model... \n'
#
##Interpolate Vs rays through this model:
#interpolation_type='linear'
#
## Run for s-waves through Qp:
#raytypes=array([1])
#Qp_ray_data=runra.interp_rays(robj,materialmodel_Qp,interpolation_type,raytypes,modtype='nodes')
##Run for s-waves through Qs:
#Qs_ray_data=runra.interp_rays(robj,materialmodel_Qs,interpolation_type,raytypes,modtype='nodes')
#
#print 'Material model made.  \n adding to the residuals object...'
#
##Add to and save in a residuals object:
##values for Qp, s-waves:
#robj.add_material_values(Qp_ray_data,materialflag_Qp,1)
#print 'Added Qp for s-waves.. \n'
##Values for Qs, s-waves:
#robj.add_material_values(Qs_ray_data,materialflag_Qs,1)
#print 'Added Qs for s-waves... \n'
#
#print 'Saving new res object...'
#rfile_interp=open(rpath_interp_Qp_Qs,'w')
#pickle.dump(robj,rfile_interp)
#rfile_interp.close()
#
#print 'Saved new res object with interpolated Q values to %s' % rpath_interp_Qp_Qs



########################################################
##############  Then compute metrics   #################
########################################################

print 'Reading residuals file...'
## Read in the file created above:
qfile = open(rpath_interp_Qp_Qs,'r')
qobj = pickle.load(qfile)
qfile.close()
print 'Read residuals file...'

print 'opening Qp... \n'
qpfile = open(hauksson_grid_path_Qp,'r')
qp_obj = pickle.load(qpfile)
qpfile.close()
print 'opened Qp model. \n'

print 'opening Qs... \n'
qsfile = open(hauksson_grid_path_Qs,'r')
qs_obj = pickle.load(qsfile)
qsfile.close()
print 'opened Qs model. \n'

## Compute metrics...
#  First for the Qp model (s-wave rays):
Qp_ray_data = qobj.rayval_s_qp

ind_s_Qp_path=ra.compute_pathintegral(Qp_ray_data,qp_obj,0)
print 'Computed path integral, Qp'
#normalization:
ind_s_Qp_normpath=ra.compute_pathintegral(Qp_ray_data,qp_obj,1)
print 'Computed normalized path integral, Qp'
#gradient:
ind_s_Qp_gradpath=ra.compute_devpathintegral(Qp_ray_data,qp_obj,0)
print 'Computed gradient integral, Qp'

# Then for Qs model (s-wave rays):
Qs_ray_data = qobj.rayval_s_qs

ind_s_Qs_path=ra.compute_pathintegral(Qs_ray_data,qs_obj,0)
print 'Computed path integral, Qs'
#normalization:
ind_s_Qs_normpath=ra.compute_pathintegral(Qs_ray_data,qs_obj,1)
print 'Computed normalized path integral, Qs'
#gradient:
ind_s_Qs_gradpath=ra.compute_devpathintegral(Qs_ray_data,qs_obj,0)
print 'Computed gradient integral, Qs'


## Then, add these to the residual object:

#indextype:          Flag for type of index: 
#                        0=path integral,1=normalized path integral, 
#                        2=gradient path integral
#ray_type:           0=P-wave, 1=S-wave
#value_flag:         Type of material model.  
#                        0=Vp, 1=Vs, 2=Vp/Vs, 3=Qp, 4=Qs, 5=Qp/Qs

# Here:
# index type =0 path, 1 norm path, 2 grad path, 
# ray type is p-0, s=1,
# material type is Qp=3, Qs=4

# First add s-wave rays through Qp, each index:
qobj.add_indices(ind_s_Qp_path,0,1,3)
#indextype =1,
qobj.add_indices(ind_s_Qp_normpath,1,1,3)
#indextype=2:
qobj.add_indices(ind_s_Qp_gradpath,2,1,3)
print 'Added S-wave metrics through Qp to object'

# Then add s-wave rays through Qs, each index:
qobj.add_indices(ind_s_Qs_path,0,1,4)
#indextype =1,
qobj.add_indices(ind_s_Qs_normpath,1,1,4)
#indextype=2:
qobj.add_indices(ind_s_Qs_gradpath,2,1,4)
print 'Added S-wave metrics through Qs to object'

print 'Saving object...'
## Save to residual object:
qfile_indices = open(rpath_indices_Qp_Qs,'w')
pickle.dump(qobj,qfile_indices)
qfile_indices.close()
print 'Qp and Qs metrics for S-waves saved into file %s' % rpath_indices_Qp_Qs
