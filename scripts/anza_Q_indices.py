#############  Script to compute Q metrics  ###########
## VJS 2/2018

import dread
import cdefs as cdf
import cPickle as pickle
from numpy import meshgrid,zeros,where,array
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

what_home=0

if what_home==0:
    #Desktop:
    HOME='/katmai'
elif what_home==1:
    #Mac:
    HOME='/Users/vsahakian'


## Path to Qp and Qs
hauksson_model_path_Qp = HOME + '/anza/data/vm/Hauksson/Qp.pckl'
hauksson_model_path_Qs = HOME + '/anza/data/vm/Hauksson/Qs.pckl'

#Path to residuals object:
rpath = HOME+'/anza/models/residuals/mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_raydat.pckl'

#Material flag - what kind is it?  0=Vp, 1=Vs, 2=Vp/Vs, 3=Qp, 4=Qs, 5=Qp/Qs
materialflag_Qp = 3
materialflag_Qs = 4

# Material model name
materialmodelname_Qp_Qs = 'HaukQp_Qs'

#Residual Object paths:
#interpolated ray values
#rpath_interp_Qp_Qs = rpath.split('.pckl')[0]+'_interp_'+materialmodelname_Qp_Qs+'.pckl'
rpath_interp_Qp_Qs = '/home/vsahakian/Desktop/' + rpath.split('.pckl')[0].split('residuals/')[1].split('/')[1] +'_interp_'+materialmodelname_Qp_Qs+'.pckl'

#interpolated ray values and indices
#rpath_indices_Qp_Qs = rpath.split('.pckl')[0]+'_interp_indices_'+materialmodelname_Qp_Qs+'.pckl'
tmp_rpath_indices_Qp_Qs = '/home/vsahakian/Desktop/' + rpath.split('.pckl')[0].split('residuals/')[1].split('/')[1] +'_interp_indices_'+materialmodelname_Qp_Qs+'.pckl'


#######################################
##############   #################
#######################################

##Interpolate rays:
#Read in residuals object
print 'opening residuals object... \n'
rfile=open(rpath,'r')
robj=pickle.load(rfile)
rfile.close()
print 'residuals object opened \n'

print 'opening Qp... \n'
qpfile = open(hauksson_model_path_Qp,'r')
materialmodel_Qp = pickle.load(qpfile)
qpfile.close()
print 'opened Qp model. \n'

print 'opening Qs... \n'
qsfile = open(hauksson_model_path_Qs,'r')
materialmodel_Qs = pickle.load(qsfile)
qsfile.close()
print 'opened Qs model. \n'

print 'Interpolating model... \n'

#Interpolate Vs rays through this model:
interpolation_type='linear'

# Run for s-waves through Qp:
raytypes=array([1])
Qp_ray_data=runra.interp_rays(robj,materialmodel_Qp,interpolation_type,raytypes)
#Run for s-waves through Qs:
Qs_ray_data=runra.interp_rays(robj,materialmodel_Qs,interpolation_type,raytypes)

print 'Material model made.  \n adding to the residuals object...'

#Add to and save in a residuals object:
#values for Qp, s-waves:
robj.add_material_values(Qp_ray_data,materialflag_Qp,1)
print 'Added Qp for s-waves.. \n'
#Values for Qs, s-waves:
robj.add_material_values(Qs_ray_data,materialflag_Qs,1)
print 'Added Qs for s-waves... \n'

print 'Saving new res object...'
rfile_interp=open(rpath_interp_Qp_Qs,'w')
pickle.dump(robj,rfile_interp)
rfile_interp.close()

print 'Saved new res object with interpolated Q values to %s' % rpath_interp_Qp_Qs
