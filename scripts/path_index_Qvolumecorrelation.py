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

what_home=0

if what_home==0:
    #Desktop:
    HOME='/katmai'
elif what_home==1:
    #Mac:
    HOME='/Users/vsahakian'

#############################
########### Paths ###########
#############################

home=HOME+'/anza/models/residuals/'
run_name = 'mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30/'
residualpath=home + run_name + 'mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_raydat.pckl'
materialpath_Qs=HOME+'/anza/data/vm/Hauksson/Qs.pckl'
materialpath_Qp=HOME+'/anza/data/vm/Hauksson/Qp.pckl'

#gobjpath=home+run_name+'/'+'mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_mean_grid.pckl'
gobjpath_mean_Qs = '/home/vsahakian/Desktop/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_meanQs_grid.pckl'
gobjpath_count_Qs = '/home/vsahakian/Desktop/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_countQs_grid.pckl'
gobjpath_median_Qs = '/home/vsahakian/Desktop/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_medianQs_grid.pckl'

gobjpath_mean_Qp = '/home/vsahakian/Desktop/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_meanQp_grid.pckl'
gobjpath_count_Qp = '/home/vsahakian/Desktop/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_countQp_grid.pckl'
gobjpath_median_Qp = '/home/vsahakian/Desktop/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_medianQp_grid.pckl'

# Set raytpe to Vs, whoich is 1:
raytype = 1

####################################################################

###############################################
########### Step 1 - Read in models ###########
###############################################
#print 'Reading in material model'
## Material model:
#mfile = open(materialpath_Qs,'r')
#mobj_Qs = pickle.load(mfile)
#mfile.close()

print 'Reading in material model'
# Material model:
mfile = open(materialpath_Qp,'r')
mobj_Qp = pickle.load(mfile)
mfile.close()

print 'Reading in residual model'
# Residuals model:
rfile = open(residualpath,'r')
robj = pickle.load(rfile)
rfile.close()

def bin_median_mean_count(robj,mobj,raytype,gobjpath_mean,gobjpath_median,gobjpath_count):
    '''
    Put the remaining steps into a function, and run separately for Qp and Qs...
    INput:
        robj:               Residuals object 
        mobj:               Material object
        raytype:            Number with type of ray. 0 = p-waves, 1=s-waves.
        gobjpath_mean:      String with path for mean
        gobjpath_median:    String with path for median
        gobjpath_count:     String with path for counts
    '''

    ###############################################
    ########### Step 2 - Get bin edges  ###########
    ###############################################
    
    # Do this all for Qs since the edges are the same as Qp
    
    print 'Getting bin edges'
    
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
    
    print 'Gridding for mean...'
    statistic = 'mean'
    gridded_obj_mean=ra.grid_path_term(robj,all_binedges,raytype,statistic,rpath_type='object',subset=ind_subset)
    
    print 'Saving mean grid...'
    # Save grid to file:
    gfile = open(gobjpath_mean,'w')
    pickle.dump(gridded_obj_mean,gfile)
    gfile.close()
    
    ####################
    print 'Gridding for counts...'
    # Also get counts:
    statistic = 'count'
    gridded_obj_count = ra.grid_path_term(robj,all_binedges,raytype,statistic,rpath_type='object',subset=ind_subset)
    
    print 'Saving count grid...'
    # Save grid to file:
    gfile = open(gobjpath_count,'w')
    pickle.dump(gridded_obj_count,gfile)
    gfile.close()
    
    ####################
    print 'Gridding for median...'
    # Also get median:
    statistic = 'median'
    gridded_obj_median = ra.grid_path_term(robj,all_binedges,raytype,statistic,rpath_type='object',subset=ind_subset)
    
    print 'Saving median grid...'
    # Save grid to file:
    gfile = open(gobjpath_median,'w')
    pickle.dump(gridded_obj_median,gfile)
    gfile.close()

###############################################
###########       Run for Qs        ###########
###############################################

#bin_median_mean_count(robj,mobj_Qs,raytype,gobjpath_mean_Qs,gobjpath_median_Qs,gobjpath_count_Qs)

bin_median_mean_count(robj,mobj_Qp,raytype,gobjpath_mean_Qp,gobjpath_median_Qp,gobjpath_count_Qp)
