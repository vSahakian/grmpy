###############################
##       Compute CI's        ## 
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

# Residuals object path:
rpath = HOME + '/anza/models/residuals/mixedregr_v3anza2013_pga_5coeff_a4_-1.72_a5_-0.01_Mc_8.5_res4_noVs30/mixedcoeff_v3anza2013_pga_noVs30_5coeff_pga__ncoeff5_Mc_8.5_VR_99.4_a4_-1.72_a5_-0.01_robj.pckl'

# Directory to save arrays:
data_dir = '/Users/vsahakian/anza/models/statistics/mixedregr_v3anza2013_pga_5coeff_a4_-1.72_a5_-0.01_Mc_8.5_res4_noVs30/'

# Directory to store figures:
fig_dir = fig_dir = '/Users/vsahakian/anza/models/statistics/mixedregr_v3anza2013_pga_5coeff_a4_-1.72_a5_-0.01_Mc_8.5_res4_noVs30/figs/'

# Basename:
bname = 'CI_1'


############
############

rfile = open(rpath,'r')
robj = pickle.load(rfile)
rfile.close()

## Get info per station:
sta_all,stnum_all,evnum_all,evlon_all,evlat_all,evdepth_all,rrup_all,path_terms_all = dbstats.extract_stationinfo(robj)

# Path std:
sigma_path = robj.path_std

# Make large CI and dPvar arrays, append all to these for later:
CI_unsort = np.array([])
dPvar_unsort = np.array([])

# For each station, compute the CI and path semivariance:
for stationi in range(len(sta_all)):
    # Compute CI...
    stai = sta_all[stationi]
    stnumi = stnum_all[stationi]
    evnumi = evnum_all[stationi]
    evloni = evlon_all[stationi]
    evlati = evlat_all[stationi]
    evdepthi = evdepth_all[stationi]
    rrupi = rrup_all[stationi]
    path_termsi = path_terms_all[stationi]
    
    print 'running station %s with %i events' % (stai,len(evnumi))
    
    # Get closeness index and dPvar for this station:
    CI_i,dPvar_i = dbstats.compute_CI_variance_perstation(stai,stnumi,evnumi,evloni,evlati,evdepthi,rrupi,path_termsi,sigma_path)
    
    # Add to the larger unsorted arrays:
    CI_unsort = np.r_[CI_unsort,CI_i]
    dPvar_unsort = np.r_[dPvar_unsort,dPvar_i]
    
    
# At the end, sort them:
CI_sort_ind = np.argsort(CI_unsort)

CI = CI_unsort[CI_sort_ind]
dPvar = dPvar_unsort[CI_sort_ind]

# Remove unnecessary arrays from memory:
del CI_unsort
del dPvar_unsort
del CI_sort_ind

# Save important ones:
datpath_ci = data_dir + bname + '_CI.pckl'
datfile = open(datpath_ci,'w')
pickle.dump(CI,datfile)
datfile.close()

datpath_dpvar = data_dir + bname + '_dPvar.pckl'
datfile = open(datpath_dpvar,'w')
pickle.dump(dPvar,datfile)
datfile.close()


# And then bin them....because this will be huge...
bins = np.linspace(0, max(CI), max(CI)*5)
digitized = np.digitize(CI, bins)
dPvar_means = [dPvar[digitized == i].mean() for i in range(1, len(bins))]

dPvar_std = [dPvar[digitized == i].std() for i in range(1,len(bins))]

    
## Plot - non nan's:
#keep_i = np.where(np.isnan(dPvar_means)==False)[0]
#x = bins[keep_i]
#y = np.array(dPvar_means)[keep_i]**2
#
#plt.figure()
#plt.scatter(x,y,facecolors='none',edgecolors='k')

plt.figure()
plt.scatter(bins[0:len(dPvar_means)],np.array(dPvar_means)**2,facecolors='none',edgecolors='k')