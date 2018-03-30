## Get transmission coefficients.
## VJS 3/2018

import numpy as np
import cPickle as pickle
import cdefs as cdf
import run_res_analysis as runra
import raytracing as rt

####################################################################
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

############################################################################
###########################      Paths       ###############################
############################################################################

# Paths to velocity model pckl, and output gridded and gradient pckls:
velmodelpath = HOME + '/anza/data/pckl/FangVs.pckl'
velmodel_gridded_path = HOME + '/anza/data/pckl/FangVs_transmissionregridded.pckl'
velmodel_gradient_path = HOME + '/anza/data/pckl/FangVs_gradient.pckl'

# Path to outptut Fang Vp and gridded Vp that will be made:
velmodelpath_vp = HOME + '/anza/data/pckl/FangVp.pckl'
velmodel_gridded_path_vp = HOME + '/anza/data/pckl/FangVp_gradient.pckl'

# Path to output density model:
denmodelpath = HOME + '/anza/data/pckl/FangDensity_brocher.pckl'

# Path to coords path and material ascii data for fang Vp
coordspath_vp = HOME+'/anza/data/vm/Fang2016/coords.txt'
materialasciipath_vp = HOME+'/anza/data/vm/Fang2016/Vp.dat'


residualpath = HOME + '/anza/models/residuals/mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_raydat.pckl'

new_z = np.array([-2., -1., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20.])

############################################################################
###################    Step 1: Read in Models    ###########################
############################################################################

# Open residuals object
rfile = open(residualpath,'r')
robj = pickle.load(rfile)
rfile.close()

# Open Fang model:
mfile = open(velmodelpath,'r')
mobj = pickle.load(mfile)
mfile.close()


############################################################################
###################  Step 2: Make Fang Vp mod    ###########################
############################################################################


#Make material model:
#Convert the longitudes from positive west to negative:
lonconvert=2
FangVp_model = runra.make_material_object(coordspath_vp,materialasciipath_vp,velmodelpath_vp,lonconvert)


############################################################################
###################  Step 3: Regrid vel model    ###########################
############################################################################

# Regrid Fang Vs and Vp according to the new depths listed in paths section
FangVs_regridded = rt.regrid_z_materialmodel(mobj,new_z)
FangVp_regridded = rt.regrid_z_materialmodel(FangVp_model,new_z)

# Save Vp to file
new_gridfile = open(velmodel_gridded_path,'w')
pickle.dump(FangVs_regridded,new_gridfile)
new_gridfile.close()

# Save Vs to file
new_gridfile = open(velmodel_gridded_path_vp,'w')
pickle.dump(FangVp_regridded,new_gridfile)
new_gridfile.close()

# Convert Vp to density using Brocher:
Density = rt.convert_velocity2density(FangVp_regridded)

## Save to file:
denfile = open(denmodelpath,'w')
pickle.dump(Density,denfile)
denfile.close()

############################################################################
###################  Step 4: Compute gradient    ###########################
############################################################################

# Need dx, dy, and dz:
dz = np.diff(FangVs_regridded.z)
dy = np.diff(FangVs_regridded.y)
dx = np.diff(FangVs_regridded.x)

# Also need the number of points:
new_nx = mobj.nx
new_ny = mobj.ny
new_nz = len(new_z)

# Set up gradient - give the sample intervals in the order of the array's shape
#   i.e., FangVs_gradient.shape = (nz,ny,nx), therefore give dz, dy, dx
FangVs_gradient_array = np.gradient(FangVs_regridded.materials,dz[0],dy[0],dx[0])

## Make a material model of it, write to file:
FangVs_gradient = cdf.material_model(mobj.x,mobj.y,new_z,new_nx,new_ny,new_nz,FangVs_gradient_array)

gfile = open(velmodel_gradient_path,'w')
pickle.dump(FangVs_gradient,gfile)
gfile.close()

############################################################################
###############  Step 5: Compute transmission per ray  #####################
############################################################################

# Loop over each ray, and compute for each ray.

## Initiate array with all values of transmission:
totaltransmission = np.array([])

#for i_ray in range(len(robj.vs_lon)):
for i_ray in range(2000):

    ## Pritn the ray number the loop is on:
    if i_ray % 1000 == 0:
        print 'Running for ray number %i' % i_ray

    ## Get the vectors/norms at every point along this ray:
    i_ray_vectors,i_ray_norms = rt.compute_raypt_vector(robj.vs_lon[i_ray],robj.vs_lat[i_ray],robj.vs_depth[i_ray])
    
    ## Get the closest index in the gradient model:
    i_ray_closestindex = rt.find_nearest_point(robj.vs_lon[i_ray],robj.vs_lat[i_ray],robj.vs_depth[i_ray],FangVs_gradient)
    
    ## Get the angle of incidence:
    i_ray_incidence_angles = rt.compute_angle_of_incidence(i_ray_vectors,i_ray_norms,i_ray_closestindex,FangVs_gradient)

    ## Get the difference across interface:
    i_ray_veldifference,i_ray_dendifference = rt.obtain_interface_properties(i_ray_closestindex,FangVs_regridded,Density)
    
    ## Compute the transmission coefficient:
    i_ray_transmission_coeff = rt.compute_transmission_coefficient(i_ray_incidence_angles,i_ray_veldifference,i_ray_dendifference)
    
    ## Compute the total transmission for this ray:
    i_ray_totaltransmission = rt.transmission_summation(i_ray_transmission_coeff)
    
    ## Append to overall array:
    totaltransmission = np.r_[totaltransmission,i_ray_totaltransmission]

############################################################################
###################  Step 5: Compute angle of inc    #######################
############################################################################

