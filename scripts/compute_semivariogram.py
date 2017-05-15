###############################
##   Compute semivariogram   ## 
###############################

# VJS 5/2017

import numpy as np
import dbstatistics as dbstats
import cPickle as pickle
import matplotlib.pyplot as plt

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

############
############

# Path to residuals object:
rpath = HOME + '/anza/models/residuals/mixedregr_v3anza2013_pga_5coeff_a4_-1.72_a5_-0.01_Mc_8.5_res4_noVs30/mixedcoeff_v3anza2013_pga_noVs30_5coeff_pga__ncoeff5_Mc_8.5_VR_99.4_a4_-1.72_a5_-0.01_robj.pckl'

# Lag distance information:
lagdistance_x = np.array([0,2,4,6,8,10,14,18,22,26,30,35,40,50,60,70,80])

# Rrup semivariogram bin information:
semivariogram_bins = np.array([15,17.5,20,22.5,25,27.5,30,32.5,35])

# Number of columns for plot:
ncols = 2

# Axis limits for plot
axlims = [[-2,70],[0,0.40]]

# Output figures directory:
fig_dir = '/Users/vsahakian/anza/models/statistics/mixedregr_v3anza2013_pga_5coeff_a4_-1.72_a5_-0.01_Mc_8.5_res4_noVs30/figs/'

# Base figure name:
bfigname = 'semivariogram1'

############
############

# Open the residuals object:
rfile = open(rpath,'r')
robj = pickle.load(rfile)
rfile.close()

####
## Bin the data according to the semivariogram bins:
evnum_all,sta_all,stnum_all,stlon_all,stlat_all,path_terms_all,rrup_all = dbstats.semivar_rbin(robj,semivariogram_bins)

###
## Loop over the bins, and per bin, compute the arrays for semivariograms.

# Set the empty list of x and y's to use later, to append x and y for each bin:
lagdistance_x_all = []
semivariance_y_all = []

# Loop:
for bin_i in range(len(semivariogram_bins)-1):
    evnum_i = evnum_all[bin_i]
    sta_i = sta_all[bin_i]
    stnum_i = stnum_all[bin_i]
    stlon_i = stlon_all[bin_i]
    stlat_i = stlat_all[bin_i]
    path_terms_i = path_terms_all[bin_i]
    rrup_i = rrup_all[bin_i]
    
    # Compute semivariance for this bin:
    lagdistance_x_i, semivariance_y_i = dbstats.semivar_variables(evnum_i,sta_i,stnum_i,stlon_i,stlat_i,path_terms_i,rrup_i,lagdistance_x)
    
    # Append to final bin:
    lagdistance_x_all.append(lagdistance_x_i)
    semivariance_y_all.append(semivariance_y_i)

# Plot them:
semivar_fig_all = dbstats.subplot_all_semivariograms(lagdistance_x_all,semivariance_y_all,ncols,axlims,semivariogram_bins)

# Save figure:
semivar_fig_all.savefig(fig_dir + 'pdf/' + bfigname + '.pdf')
semivar_fig_all.savefig(fig_dir + bfigname + '.png')
