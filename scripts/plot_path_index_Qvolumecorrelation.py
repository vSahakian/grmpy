#### Compre gridded path terms with a material model
## VJS 12/2017

import dread as dr
import cPickle as pickle
import res_analysis as ra
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr


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

#############################
########### Paths ###########
#############################

home=HOME+'/anza/models/residuals/'
run_name = 'mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30/'

residualpath = home + run_name + 'mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_raydat.pckl'
materialpath_Qs = HOME+'/anza/data/vm/Hauksson/Qs.pckl'

#gobjpath=home+run_name+'/'+'mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_mean_grid.pckl'

#gobjpath_mean = '/home/vsahakian/Desktop/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_mean_grid.pckl'
#gobjpath_count = '/home/vsahakian/Desktop/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_count_grid.pckl'
#gobjpath_median = '/home/vsahakian/Desktop/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_median_grid.pckl'

gobjpath_mean_Qs = home + run_name + 'mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_meanQs_grid.pckl'
gobjpath_count_Qs = home + run_name + 'mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_countQs_grid.pckl'
gobjpath_median_Qs = home + run_name + 'mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_medianQs_grid.pckl'

png_dir = home + run_name + 'figs/'
pdf_dir = home + run_name + 'figs/pdfs/'

statistic = 'mean'

# Set raytpe to Qs, whoich is 1:
raytype_Qs = 3

####################################################################

###############################################
###########  Step 1 - Read in grids ###########
###############################################
print 'Reading in material model'
mfile = open(materialpath_Qs,'r')
mobj = pickle.load(mfile)
mfile.close()

####################
# Gridded path mean:
print 'Reading in mean model'
gfile = open(gobjpath_mean_Qs,'r')
gridded_obj_mean_Qs = pickle.load(gfile)
gfile.close()

####################
# Gridded path counts:
print 'Reading in count model'
gfile = open(gobjpath_count_Qs,'r')
gridded_obj_count_Qs = pickle.load(gfile)
gfile.close()

####################
# Gridded path median:
print 'Reading in median model'
gfile = open(gobjpath_median_Qs,'r')
gridded_obj_median_Qs = pickle.load(gfile)
gfile.close()



###############################################
###########     Step 2 - Reshape    ###########
###############################################

print 'Reshaping...'
# First reshape the statistics:
# Mean
print 'Reshaping mean statistic'
# For every depth slice, convert the x by y to a vector; then concatenate all at the end
for i_mean_depth in range(np.shape(gridded_obj_mean_Qs.statistic)[2]):
    i_mean_array_length = np.shape(gridded_obj_mean_Qs.statistic)[0] * np.shape(gridded_obj_mean_Qs.statistic)[1]
    i_mean_slice = np.reshape(gridded_obj_mean_Qs.statistic[:,:,i_mean_depth],[1,i_mean_array_length])
    
    # If it's the first slice, make the mean reshape equal to it:
    if i_mean_depth == 0:
        mean_reshape_Qs = i_mean_slice
    else:
        mean_reshape_Qs = np.append(mean_reshape_Qs,i_mean_slice)

# Median
print 'Reshaping median statistic'
# For every depth slice, convert the x by y to a vector; then concatenate all at the end
for i_median_depth in range(np.shape(gridded_obj_median_Qs.statistic)[2]):
    i_median_array_length = np.shape(gridded_obj_median_Qs.statistic)[0] * np.shape(gridded_obj_median_Qs.statistic)[1]
    i_median_slice = np.reshape(gridded_obj_median_Qs.statistic[:,:,i_median_depth],[1,i_median_array_length])
    
    # If it's the first slice, make the mean reshape equal to it:
    if i_median_depth == 0:
        median_reshape_Qs = i_median_slice
    else:
        median_reshape_Qs = np.append(median_reshape_Qs,i_median_slice)

# Count
print 'Reshaping count statistic'
# For every depth slice, convert the x by y to a vector; then concatenate all at the end
for i_count_depth in range(np.shape(gridded_obj_count_Qs.statistic)[2]):
    i_count_array_length = np.shape(gridded_obj_count_Qs.statistic)[0] * np.shape(gridded_obj_count_Qs.statistic)[1]
    i_count_slice = np.reshape(gridded_obj_count_Qs.statistic[:,:,i_count_depth],[1,i_count_array_length])
    
    # If it's the first slice, make the mean reshape equal to it:
    if i_count_depth == 0:
        count_reshape_Qs = i_count_slice
    else:
        count_reshape_Qs = np.append(count_reshape_Qs,i_count_slice)

## Now reshape the velocity model 
#  For every depth slice, convert the x by y to a vector; then concatenate all at the end.
for i_depth in range(np.shape(mobj.materials)[0]):
    i_array_length = np.shape(mobj.materials)[1] * np.shape(mobj.materials)[2]
    i_slice = np.reshape(mobj.materials[i_depth,:,:], [1,i_array_length])

    # If it's the first slice, make the velocity model equal to it:
    if i_depth == 0:
        Qmodel_reshape = i_slice
    else:
        Qmodel_reshape = np.append(Qmodel_reshape,i_slice)

# Also use np.sparse

###############################################
########### Step 5 - Correlate bins  ##########
###############################################
# First, find where there are NaN's:
mean_isnan_ind = np.where(np.isnan(mean_reshape_Qs) == True)[0]
median_isnan_ind = np.where(np.isnan(median_reshape_Qs) == True)[0]

# Find where counts is 0:
count_0_ind = np.where(count_reshape_Qs == 0)[0]

# If they are the same:
if len(np.where((mean_isnan_ind - count_0_ind) !=0)[0]) == 0:
    print 'Mean NaNs and count 0s are the same cells'
if len(np.where((median_isnan_ind - count_0_ind) !=0)[0]) == 0:
    print 'Median NaNs and count 0s are the same cells'

## 
# Now, find where they are NOT NaN:
mean_Qs_nonan_ind = np.where(np.isnan(mean_reshape_Qs) == False)[0]
median_Qs_nonan_ind = np.where(np.isnan(median_reshape_Qs) == False)[0]

# Get only the non-Nan values:
mean_Qs_reshape_noNan = mean_reshape_Qs[mean_Qs_nonan_ind]
median_Qs_reshape_noNan = median_reshape_Qs[median_Qs_nonan_ind]

# Now also get the Qmodel at these indices:
Qs_model_reshape_noNan = Qmodel_reshape[mean_Qs_nonan_ind]

## Get the pearsonr:
mean_Qs_r,mean_Qs_p = pearsonr(mean_Qs_reshape_noNan,Qs_model_reshape_noNan)
median_Qs_r,median_Qs_p = pearsonr(median_Qs_reshape_noNan,Qs_model_reshape_noNan)

###############################################
###########     Step 6 - Plots      ###########
###############################################

# First make a scatter plot of velocity and mean path term
scatterfig = plt.figure(figsize=(12,9))

# Scatter:
plt.scatter(Qs_model_reshape_noNan,mean_Qs_reshape_noNan,marker='o',color='#2f2772',label='Mean path term - r = %.2f, p = %.1f' % (mean_Qs_r,mean_Qs_p))
plt.scatter(Qs_model_reshape_noNan,median_Qs_reshape_noNan,marker='^',color='#30705b',label='Median path term - r = %.2f, p = %.1f' % (median_Qs_r,median_Qs_p))

plt.ylim([-2.5,2])
plt.xlim([0,1200])

plt.xlabel('Cell Qs (km/s)')
plt.ylabel('Path Term (ln residual)')
plt.title('Gridded path term vs. Qs')

plt.legend(loc=1)

# Save figure:
figname = 'griddedpath_vs_Qs'
scatterfig.savefig(png_dir + figname + '.png')
scatterfig.savefig(pdf_dir + figname + '.pdf')


###########################################################
## Make slice plots through vel model and path model
slicefig,(velax,pathax) = plt.subplots(nrows=1,ncols=2,figsize=(15,9))

# First plot the velocity model slice:
mobj.plot_xslice(velax,-116,'magma',[0,1150],'Latitude','Depth (km)','HaukssonQs',0.1,[[32.5,34],[-2,22]])

# Then the path axis:
gridded_obj_mean_Qs.plot_slice(pathax,-116,'lon',0.1,'magma',[-2,2],[[32.5,34],[0,22]],binnodes=np.array([mobj.x,mobj.y,mobj.z]))