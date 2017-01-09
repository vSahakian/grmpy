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

what_home=0

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

# Model 1 - Inversion of Anza data
#   Residuals object:
tradpath=home+'v2anza2013_Mc8.5_pgrid_30sta_res4/anza2013_Mc8.5_pgrid_30sta_robj.pckl'

# Model 2 - Mixed Effects inversion of Anza data
mixedpath=home+'mixedregr_v2anza2013_Mc_8.5_30sta/mixedregr_v2anza2013_30sta_Mc8.5_VR_99.9_robj.pckl'

# Comparison directory name:
comp_dir = 'compare_v2anza2013_30sta'

# Plotting colormap:
mymap='gnuplot'

evaxlim = [[-6,6],[-6,6]]
paxlim = [[-6,6],[-6,6]]
staxlim = [[-6,6],[-6,6]]

##############################################################################

# Initiate directory:
ra.compare_init(home,comp_dir)

# Plot:
ra.compare_mixed_traditional(home,comp_dir,tradpath,mixedpath,mymap,evaxlim,staxlim,paxlim)


