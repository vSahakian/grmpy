#### Compre gridded path terms with a material model
## VJS 12/2017

import dread as dr
import cPickle as pickle
import res_analysis as ra
import numpy as np
import matplotlib.pyplot as plt


####################################################################
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

#############################
########### Paths ###########
#############################

home=HOME+'/anza/models/residuals/'
run_name = 'mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30/'
residualpath=home + run_name + 'mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_raydat.pckl'
materialpath=HOME+'/anza/data/pckl/FangVs.pckl'
gobjpath=home+run_name+'/'+'mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_mean_grid.pckl'

statistic = 'mean'

# Set raytpe to Vs, whoich is 1:
raytype = 1

####################################################################

###############################################
########### Step 1 - Read in models ###########
###############################################

# Material model:
mfile = open(materialpath,'r')
mobj = pickle.load(mfile)
mfile.close()

# Residuals model:
rfile = open(residualpath,'r')
robj = pickle.load(rfile)
rfile.close()


###############################################
########### Step 2 - Get bin edges  ###########
###############################################

# Material object has nodes - that means that there will be len(nodes) + 2 edges:
x_numedges = len(mobj.x) + 1
y_numedges = len(mobj.y) + 1
z_numedges = len(mobj.z) + 1

# Get the differences between nodes:
x_diff = np.diff(mobj.x)
y_diff = np.diff(mobj.y)
z_diff = np.diff(mobj.z)

# Initiate empty bin edges:
x_binedges = np.zeros((x_numedges))
y_binedges = np.zeros((y_numedges))
z_binedges = np.zeros((z_numedges))


# Fill first of each:
x_binedges[0] = mobj.x[0] - (x_diff[0]/2.)
y_binedges[0] = mobj.y[0] - (y_diff[0]/2.)
z_binedges[0] = mobj.z[0] - (z_diff[0]/2.)

#
for i_node in range(len(mobj.x)):
    if i_node < len(mobj.x) - 1:
        x_binedges[i_node + 1] = mobj.x[i_node] + (x_diff[i_node]/2.)
    else:
        x_binedges[i_node + 1] = mobj.x[i_node] + (x_diff[i_node - 1]/2.)
    
for i_node in range(len(mobj.y)):
    if i_node < len(mobj.y) - 1:
        y_binedges[i_node + 1] = mobj.y[i_node] + (y_diff[i_node]/2.)
    else:
        y_binedges[i_node + 1] = mobj.y[i_node] + (y_diff[i_node - 1]/2.)

for i_node in range(len(mobj.z)):
    if i_node < len(mobj.z) - 1:
        z_binedges[i_node + 1] = mobj.z[i_node] + (z_diff[i_node]/2.)
    else:
        z_binedges[i_node + 1] = mobj.z[i_node] + (z_diff[i_node - 1]/2.)

# Make binedges arrays:
all_binedges = [x_binedges,y_binedges,z_binedges]


###############################################
########### Step 3 - Compute grid   ###########
###############################################

# Use all rays:
ind_subset = 'all'

gridded_obj=ra.grid_path_term(robj,all_binedges,raytype,statistic,rpath_type='object',subset=ind_subset)


###############################################
###########     Step 4 - Reshape    ###########
###############################################



###############################################
########### Step 5 - Correlate bins  ##########
###############################################


###############################################
###########     Step 6 - Plots      ###########
###############################################



