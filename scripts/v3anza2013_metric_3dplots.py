## Scatter plot of residuals, metric, and r in 3d
## Test to see what surface to fit to it to get coefficients...
# VJS 10/2017

import cPickle as pickle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

what_home=1

if what_home==0:
    #Desktop:
    HOME='/media/vsahakian/katmai'
elif what_home==1:
    #Mac:
    HOME='/Users/vsahakian'

##### Paths #####

# Path to residuals object:
rpath=HOME+'/anza/models/residuals/mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_raydat_interp_indices_FangVs.pckl'


##################

# Open pickle:
rfile = open(rpath,'r')
robj = pickle.load(rfile)
rfile.close()

##################

# Start 3d plot:

# without normalization...

fig = plt.figure(figsize=(40,30))
ax3d = fig.add_subplot(1,1,1,axisbg="1.0")
ax3d = fig.gca(projection='3d')

ax3d.scatter(robj.ind_s_vs_gradpathint,robj.r,robj.path_terms,s=4,depthshade=True)
ax3d.set_xlabel('Gradient metric')
ax3d.set_ylabel('Rrup')
ax3d.set_zlabel('Path Residual')


# with distance normalization...

fig = plt.figure(figsize=(40,30))
ax3d = fig.add_subplot(1,1,1,axisbg="1.0")
ax3d = fig.gca(projection='3d')

ax3d.scatter(robj.ind_s_vs_gradpathint/(robj.r/max(robj.r)),robj.r,robj.path_terms,s=4,depthshade=True)
ax3d.set_xlabel('Gradient metric')
ax3d.set_ylabel('Rrup')
ax3d.set_zlabel('Path Residual')


# with distance normalization by just distance...

fig = plt.figure(figsize=(40,30))
ax3d = fig.add_subplot(1,1,1,axisbg="1.0")
ax3d = fig.gca(projection='3d')

ax3d.scatter(robj.ind_s_vs_gradpathint/robj.r,robj.r,robj.path_terms,s=4,depthshade=True)
ax3d.set_xlabel('Gradient metric')
ax3d.set_ylabel('Rrup')
ax3d.set_zlabel('Path Residual')

