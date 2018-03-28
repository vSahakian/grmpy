## Get transmission coefficients.
## VJS 3/2018

import numpy as np
import cPickle as pickle
from scipy import interpolate
import cdefs as cdf

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


velmodelpath = HOME + '/anza/data/pckl/FangVs.pckl'
velmodel_gridded_path = HOME + '/anza/data/pckl/FangVs_transmissionregridded.pckl'
velmodel_gradient_path = HOME + '/anza/data/pckl/FangVs_gradient.pckl'


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
###################  Step 2: Regrid Vel model    ###########################
############################################################################

# First get current points:
current_lon_arr = mobj.x
current_lat_arr = mobj.y
current_dep_arr = mobj.z

# Meshgrid them, to have each point:
current_dep,current_lat,current_lon = np.meshgrid(current_dep_arr,current_lat_arr,current_lon_arr,indexing='ij')

# Meshgrid the new points, with the new z array:
new_dep,new_lat,new_lon = np.meshgrid(new_z,current_lat_arr,current_lon_arr,indexing='ij')

# Get numbres of points:
current_nx = mobj.nx
current_ny = mobj.ny
current_nz = mobj.nz

new_nx = mobj.nx
new_ny = mobj.ny
new_nz = len(new_z)

# Reshape everything, by each z-slice:
lon_points = np.array([])
lat_points = np.array([])
dep_points = np.array([])
vs_values = np.array([])

lon_interp = np.array([])
lat_interp = np.array([])
dep_interp = np.array([])

for i_depth in range(current_nz):
    # Reshape the existing values of everything:
    i_lon_pts = current_lon[i_depth].ravel()
    i_lat_pts = current_lat[i_depth].ravel()
    i_dep_pts = current_dep[i_depth].ravel()
    
    i_vs_vals = mobj.materials[i_depth].ravel()
    
    # Append these to the arrays:
    lon_points = np.r_[lon_points,i_lon_pts]
    lat_points = np.r_[lat_points,i_lat_pts]
    dep_points = np.r_[dep_points,i_dep_pts]
    
    vs_values = np.r_[vs_values,i_vs_vals]
    
# And now for interp points, to get them as arrays:
for i_depth in range(new_nz):
    
    i_lon_interp = new_lon[i_depth].ravel()
    i_lat_interp = new_lat[i_depth].ravel()
    i_dep_interp = new_dep[i_depth].ravel()
    
    # Append these to the arrays:    
    lon_interp = np.r_[lon_interp,i_lon_interp]
    lat_interp = np.r_[lat_interp,i_lat_interp]
    dep_interp = np.r_[dep_interp,i_dep_interp]
    
## Re grid them to get vs_values at new points
# concatenate points:
existing_points = np.c_[lon_points,lat_points,dep_points]
# grid:
vs_interp_array = interpolate.griddata(existing_points,vs_values,(lon_interp,lat_interp,dep_interp),method='linear')

## Re-format from 1D arrays to properly sized 3D array:
vs_interp = np.zeros((new_nz,new_ny,new_nx))
for i_depth in range(new_nz):
    # Reshape to new length of x and y, for each slice:
    slice_length = new_nx * new_ny
    slice_shape = (new_ny,new_nx)
    vs_interp[i_depth] = vs_interp_array[i_depth*slice_length:(i_depth+1)*slice_length].reshape(slice_shape)
    
## Make a new velocity model object, and save to file:
FangVs_regridded = cdf.material_model(mobj.x,mobj.y,new_z,new_nx,new_ny,new_nz,vs_interp)

new_gridfile = open(velmodel_gridded_path,'w')
pickle.dump(FangVs_regridded,new_gridfile)
new_gridfile.close()

############################################################################
###################  Step 3: Compute gradient    ###########################
############################################################################

# Need dx, dy, and dz:
dz = np.diff(FangVs_regridded.z)
dy = np.diff(FangVs_regridded.y)
dx = np.diff(FangVs_regridded.x)

# Set up gradient - give the sample intervals in the order of the array's shape
#   i.e., FangVs_gradient.shape = (nz,ny,nx), therefore give dz, dy, dx
FangVs_gradient_array = np.gradient(FangVs_regridded.materials,dz[0],dy[0],dx[0])

## Make a material model of it, write to file:
FangVs_gradient = cdf.material_model(mobj.x,mobj.y,new_z,new_nx,new_ny,new_nz,FangVs_gradient_array)

gfile = open(velmodel_gradient_path,'w')
pickle.dump(FangVs_gradient,gfile)
gfile.close()

############################################################################
##################  Step 4: Get grad vec norms   ###########################
############################################################################



############################################################################
###################  Step 5: Compute angle of inc    #######################
############################################################################

