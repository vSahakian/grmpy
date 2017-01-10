# Compare event, site, and path terms between models
import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import res_analysis as ra

#Change this parameter depending on where you run:
#0=desktop
#1=mac

what_home=1

if what_home==0:
    #Desktop:
    HOME='/media/vsahakian/katmai'
    codehome='/home/vsahakian'
elif what_home==1:
    #Mac:
    HOME='/Users/vsahakian'
    codehome='/Users/vsahakian'

# Set homes:
home=HOME+'/anza/models/residuals/'

#################################################################################
############## Comparison #1 - minimum of 30 stations per event #################
#################################################################################

# Model 1 - Inversion of Anza data
#   Residuals object:
tradpath=home+'v2anza2013_Mc8.5_pgrid_30sta_res4/v2anza2013_Mc8.5_pgrid_30sta_res4_robj.pckl'

# Model 2 - Mixed Effects inversion of Anza data
mixedpath=home+'mixedregr_v2anza2013_Mc_8.5_30sta/mixedregr_v2anza2013_30sta_Mc_8.5_VR_99.9_robj.pckl'

# Comparison directory name:
comp_dir = 'compare_v2anza2013_30sta'

# Plotting colormap:
mymap='gnuplot'

evaxlim = [[-1.2,1.2],[-1.2,1.2]]
paxlim = [[-3,3],[-3,3]]
staxlim = [[-2,2],[-2,2]]

# Symbol size: symbol_size = [event_size, path_size, site_size]
symbol_size = [20,20,20]

# Color bin size for colorscale: cbins = [event,site]
cbins = [0.01, 0.01]
##############################################################################

# Initiate directory:
ra.compare_init(home,comp_dir)

# Plot:
ra.compare_mixed_traditional(home,comp_dir,tradpath,mixedpath,mymap,evaxlim,staxlim,paxlim,symbol_size,cbins)




#################################################################################
############## Comparison #2 - minimum of 32 stations per event #################
#################################################################################
# Model 1 - Inversion of Anza data
#   Residuals object:
tradpath=home+'v2anza2013_Mc8.5_pgrid_32sta_res4/v2anza2013_Mc8.5_pgrid_32sta_res4_robj.pckl'

# Model 2 - Mixed Effects inversion of Anza data
mixedpath=home+'mixedregr_v2anza2013_Mc_8.5_32sta/mixedregr_v2anza2013_32sta_Mc_8.5_VR_99.9_robj.pckl'

# Comparison directory name:
comp_dir = 'compare_v2anza2013_32sta'

# Plotting colormap:
mymap='gnuplot'

evaxlim = [[-1.2,1.2],[-1.2,1.2]]
paxlim = [[-3,3],[-3,3]]
staxlim = [[-2,2],[-2,2]]

# Symbol size: symbol_size = [event_size, path_size, site_size]
symbol_size = [20,20,20]

# Color bin size for colorscale: cbins = [event,site]
cbins = [0.01, 0.01]
##############################################################################

# Initiate directory:
ra.compare_init(home,comp_dir)

# Plot:
ra.compare_mixed_traditional(home,comp_dir,tradpath,mixedpath,mymap,evaxlim,staxlim,paxlim,symbol_size,cbins)