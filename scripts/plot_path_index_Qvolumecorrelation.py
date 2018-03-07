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
materialpath_Qp = HOME+'/anza/data/vm/Hauksson/Qp.pckl'

#gobjpath=home+run_name+'/'+'mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_mean_grid.pckl'

#gobjpath_mean = '/home/vsahakian/Desktop/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_mean_grid.pckl'
#gobjpath_count = '/home/vsahakian/Desktop/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_count_grid.pckl'
#gobjpath_median = '/home/vsahakian/Desktop/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_median_grid.pckl'

gobjpath_mean_Qs = home + run_name + 'mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_meanQs_grid.pckl'
gobjpath_count_Qs = home + run_name + 'mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_countQs_grid.pckl'
gobjpath_median_Qs = home + run_name + 'mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_medianQs_grid.pckl'

gobjpath_mean_Qp = home + run_name + 'mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_meanQp_grid.pckl'
gobjpath_count_Qp = home + run_name + 'mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_countQp_grid.pckl'
gobjpath_median_Qp = home + run_name + 'mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_medianQp_grid.pckl'

png_dir = home + run_name + 'figs/'
pdf_dir = home + run_name + 'figs/pdfs/'

statistic = 'mean'

# Set raytpe to Qs, whoich is 1:
raytype_Qs = 3

####################################################################

###############################################
###########  Step 1 - Read in grids ###########
###############################################
print 'Reading in material Qs model'
mfile = open(materialpath_Qs,'r')
mobj_Qs = pickle.load(mfile)
mfile.close()

print 'Reading in material Qp model'
mfile = open(materialpath_Qp,'r')
mobj_Qp = pickle.load(mfile)
mfile.close()

####################
# Gridded path mean:
print 'Reading in mean model, Qs'
gfile = open(gobjpath_mean_Qs,'r')
gridded_obj_mean_Qs = pickle.load(gfile)
gfile.close()

print 'Reading in mean model, Qp'
gfile = open(gobjpath_mean_Qp,'r')
gridded_obj_mean_Qp = pickle.load(gfile)
gfile.close()


####################
# Gridded path counts:
print 'Reading in count model, Qs'
gfile = open(gobjpath_count_Qs,'r')
gridded_obj_count_Qs = pickle.load(gfile)
gfile.close()

print 'Reading in count model, Qp'
gfile = open(gobjpath_count_Qp,'r')
gridded_obj_count_Qp = pickle.load(gfile)
gfile.close()

####################
# Gridded path median:
print 'Reading in median model, Qs'
gfile = open(gobjpath_median_Qs,'r')
gridded_obj_median_Qs = pickle.load(gfile)
gfile.close()

print 'Reading in median model, Qp'
gfile = open(gobjpath_median_Qp,'r')
gridded_obj_median_Qp = pickle.load(gfile)
gfile.close()


###############################################
###########     Step 2 - Reshape    ###########
###############################################

print 'Reshaping...'
# First reshape the statistics:
# Mean
print 'Reshaping mean statistic, Qs...'
# For every depth slice, convert the x by y to a vector; then concatenate all at the end
for i_mean_depth in range(np.shape(gridded_obj_mean_Qs.statistic)[2]):
    i_mean_array_length = np.shape(gridded_obj_mean_Qs.statistic)[0] * np.shape(gridded_obj_mean_Qs.statistic)[1]
    i_mean_slice = np.reshape(gridded_obj_mean_Qs.statistic[:,:,i_mean_depth],[1,i_mean_array_length])
    
    # If it's the first slice, make the mean reshape equal to it:
    if i_mean_depth == 0:
        mean_reshape_Qs = i_mean_slice
    else:
        mean_reshape_Qs = np.append(mean_reshape_Qs,i_mean_slice)
        
print 'Reshaping mean statistic, Qp'
# For every depth slice, convert the x by y to a vector; then concatenate all at the end
for i_mean_depth in range(np.shape(gridded_obj_mean_Qp.statistic)[2]):
    i_mean_array_length = np.shape(gridded_obj_mean_Qp.statistic)[0] * np.shape(gridded_obj_mean_Qp.statistic)[1]
    i_mean_slice = np.reshape(gridded_obj_mean_Qp.statistic[:,:,i_mean_depth],[1,i_mean_array_length])
    
    # If it's the first slice, make the mean reshape equal to it:
    if i_mean_depth == 0:
        mean_reshape_Qp = i_mean_slice
    else:
        mean_reshape_Qp = np.append(mean_reshape_Qp,i_mean_slice)
        
###############################################
# Median
print 'Reshaping median statistic, Qs'
# For every depth slice, convert the x by y to a vector; then concatenate all at the end
for i_median_depth in range(np.shape(gridded_obj_median_Qs.statistic)[2]):
    i_median_array_length = np.shape(gridded_obj_median_Qs.statistic)[0] * np.shape(gridded_obj_median_Qs.statistic)[1]
    i_median_slice = np.reshape(gridded_obj_median_Qs.statistic[:,:,i_median_depth],[1,i_median_array_length])
    
    # If it's the first slice, make the mean reshape equal to it:
    if i_median_depth == 0:
        median_reshape_Qs = i_median_slice
    else:
        median_reshape_Qs = np.append(median_reshape_Qs,i_median_slice)
        
print 'Reshaping median statistic, Qp'
# For every depth slice, convert the x by y to a vector; then concatenate all at the end
for i_median_depth in range(np.shape(gridded_obj_median_Qp.statistic)[2]):
    i_median_array_length = np.shape(gridded_obj_median_Qp.statistic)[0] * np.shape(gridded_obj_median_Qp.statistic)[1]
    i_median_slice = np.reshape(gridded_obj_median_Qp.statistic[:,:,i_median_depth],[1,i_median_array_length])
    
    # If it's the first slice, make the mean reshape equal to it:
    if i_median_depth == 0:
        median_reshape_Qp = i_median_slice
    else:
        median_reshape_Qp = np.append(median_reshape_Qp,i_median_slice)

###############################################
# Count
print 'Reshaping count statistic, Qs'
# For every depth slice, convert the x by y to a vector; then concatenate all at the end
for i_count_depth in range(np.shape(gridded_obj_count_Qs.statistic)[2]):
    i_count_array_length = np.shape(gridded_obj_count_Qs.statistic)[0] * np.shape(gridded_obj_count_Qs.statistic)[1]
    i_count_slice = np.reshape(gridded_obj_count_Qs.statistic[:,:,i_count_depth],[1,i_count_array_length])
    
    # If it's the first slice, make the mean reshape equal to it:
    if i_count_depth == 0:
        count_reshape_Qs = i_count_slice
    else:
        count_reshape_Qs = np.append(count_reshape_Qs,i_count_slice)
        
        
print 'Reshaping count statistic, Qp'
# For every depth slice, convert the x by y to a vector; then concatenate all at the end
for i_count_depth in range(np.shape(gridded_obj_count_Qp.statistic)[2]):
    i_count_array_length = np.shape(gridded_obj_count_Qp.statistic)[0] * np.shape(gridded_obj_count_Qp.statistic)[1]
    i_count_slice = np.reshape(gridded_obj_count_Qp.statistic[:,:,i_count_depth],[1,i_count_array_length])
    
    # If it's the first slice, make the mean reshape equal to it:
    if i_count_depth == 0:
        count_reshape_Qp = i_count_slice
    else:
        count_reshape_Qp = np.append(count_reshape_Qp,i_count_slice)
        

###############################################
## Now reshape the velocity model 
#  For every depth slice, convert the x by y to a vector; then concatenate all at the end.
print 'Reshaping the Qs Hauksson grid...'
for i_depth in range(np.shape(mobj_Qs.materials)[0]):
    i_array_length = np.shape(mobj_Qs.materials)[1] * np.shape(mobj_Qs.materials)[2]
    i_slice = np.reshape(mobj_Qs.materials[i_depth,:,:], [1,i_array_length])

    # If it's the first slice, make the velocity model equal to it:
    if i_depth == 0:
        Qsmodel_reshape = i_slice
    else:
        Qsmodel_reshape = np.append(Qsmodel_reshape,i_slice)
        
print 'Reshaping the Qp Hauksson grid...'
for i_depth in range(np.shape(mobj_Qp.materials)[0]):
    i_array_length = np.shape(mobj_Qp.materials)[1] * np.shape(mobj_Qp.materials)[2]
    i_slice = np.reshape(mobj_Qp.materials[i_depth,:,:], [1,i_array_length])

    # If it's the first slice, make the velocity model equal to it:
    if i_depth == 0:
        Qpmodel_reshape = i_slice
    else:
        Qpmodel_reshape = np.append(Qpmodel_reshape,i_slice)


###############################################
########### Step 5 - Correlate bins  ##########
###############################################
## First, find where there are NaN's:
mean_isnan_ind_Qs = np.where(np.isnan(mean_reshape_Qs) == True)[0]
median_isnan_ind_Qs = np.where(np.isnan(median_reshape_Qs) == True)[0]

mean_isnan_ind_Qp = np.where(np.isnan(mean_reshape_Qp) == True)[0]
median_isnan_ind_Qp = np.where(np.isnan(median_reshape_Qp) == True)[0]

## Find where counts is 0:
count_0_ind_Qs = np.where(count_reshape_Qs == 0)[0]
count_0_ind_Qp = np.where(count_reshape_Qp == 0)[0]

## If they are the same:
# Qs:
if len(np.where((mean_isnan_ind_Qs - count_0_ind_Qs) !=0)[0]) == 0:
    print 'Qs Mean NaNs and count 0s are the same cells'
if len(np.where((median_isnan_ind_Qs - count_0_ind_Qs) !=0)[0]) == 0:
    print 'Qs Median NaNs and count 0s are the same cells'
    
# Qp:
if len(np.where((mean_isnan_ind_Qp - count_0_ind_Qp) !=0)[0]) == 0:
    print 'Qp Mean NaNs and count 0s are the same cells'
if len(np.where((median_isnan_ind_Qp - count_0_ind_Qp) !=0)[0]) == 0:
    print 'Qp Median NaNs and count 0s are the same cells'

#### 
## Now, find where they are NOT NaN:
# Qs:
mean_Qs_nonan_ind = np.where(np.isnan(mean_reshape_Qs) == False)[0]
median_Qs_nonan_ind = np.where(np.isnan(median_reshape_Qs) == False)[0]

# Qp:
mean_Qp_nonan_ind = np.where(np.isnan(mean_reshape_Qp) == False)[0]
median_Qp_nonan_ind = np.where(np.isnan(median_reshape_Qp) == False)[0]


## Get only the non-Nan values:
# Qs:
mean_Qs_reshape_noNan = mean_reshape_Qs[mean_Qs_nonan_ind]
median_Qs_reshape_noNan = median_reshape_Qs[median_Qs_nonan_ind]

# Qp:
mean_Qp_reshape_noNan = mean_reshape_Qp[mean_Qp_nonan_ind]
median_Qp_reshape_noNan = median_reshape_Qp[median_Qp_nonan_ind]

## Now also get the Qmodel at these indices:
# Qs:
Qs_model_reshape_noNan = Qsmodel_reshape[mean_Qs_nonan_ind]

# Qp:
Qp_model_reshape_noNan = Qpmodel_reshape[mean_Qp_nonan_ind]

## Get the pearsonr:
# Qs:
mean_Qs_r,mean_Qs_p = pearsonr(mean_Qs_reshape_noNan,Qs_model_reshape_noNan)
median_Qs_r,median_Qs_p = pearsonr(median_Qs_reshape_noNan,Qs_model_reshape_noNan)

# Qp:
mean_Qp_r,mean_Qp_p = pearsonr(mean_Qp_reshape_noNan,Qp_model_reshape_noNan)
median_Qp_r,median_Qp_p = pearsonr(median_Qp_reshape_noNan,Qp_model_reshape_noNan)

###############################################
###########     Step 6 - Plots      ###########
###############################################

### Qs ####

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


### Qp ####

# First make a scatter plot of velocity and mean path term
scatterfig = plt.figure(figsize=(12,9))

# Scatter:
plt.scatter(Qp_model_reshape_noNan,mean_Qp_reshape_noNan,marker='o',color='#2f2772',label='Mean path term - r = %.2f, p = %.1f' % (mean_Qp_r,mean_Qp_p))
plt.scatter(Qp_model_reshape_noNan,median_Qp_reshape_noNan,marker='^',color='#30705b',label='Median path term - r = %.2f, p = %.1f' % (median_Qp_r,median_Qp_p))

plt.ylim([-2.5,2])
plt.xlim([0,1200])

plt.xlabel('Cell Qp (km/s)')
plt.ylabel('Path Term (ln residual)')
plt.title('Gridded path term vs. Qp')

plt.legend(loc=1)

# Save figure:
figname = 'griddedpath_vs_Qp'
scatterfig.savefig(png_dir + figname + '.png')
scatterfig.savefig(pdf_dir + figname + '.pdf')



############################################################
### Make slice plots through vel model and path model
#slicefig,(velax,pathax) = plt.subplots(nrows=1,ncols=2,figsize=(15,9))
#
## First plot the velocity model slice:
#mobj_Qs.plot_xslice(velax,-116,'magma',[0,1150],'Latitude','Depth (km)','HaukssonQs',0.1,[[32.5,34],[-2,22]])
#
## Then the path axis:
#gridded_obj_mean_Qs.plot_slice(pathax,-116,'lon',0.1,'magma',[-2,2],[[32.5,34],[0,22]],binnodes=np.array([mobj_Qs.x,mobj_Qs.y,mobj_Qs.z]))