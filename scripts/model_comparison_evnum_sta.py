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
############## Comparison #0 - minimum of 5 stations per event ##################
#################################################################################
# Model 1 - Inversion of Anza data
#   Residuals object:
tradpath_5=home+'v2anza2013_Mc8.5_pgrid_5sta_res4/v2anza2013_Mc8.5_pgrid_5sta_res4_robj.pckl'

# Model 2 - Mixed Effects inversion of Anza data
mixedpath_5=home+'mixedregr_v2anza2013_Mc_8.5_res4/mixedregr_v2anza2013_Mc_8.5_VR_99.9_robj.pckl'

# Comparison directory name:
comp_dir = 'compare_v2anza2013_5sta'

# Plotting colormap:
mymap='gnuplot'

evaxlim = [[-4.5,4.5],[-4.5,4.5]]
paxlim = [[-4.5,4.5],[-4.5,4.5]]
staxlim = [[-4.5,4.5],[-4.5,4.5]]

# Symbol size: symbol_size = [event_size, path_size, site_size]
symbol_size = [15,15,15]

# Color bin size for colorscale: cbins = [event,site]
cbins = [0.01, 10]
##############################################################################

# Initiate directory:
ra.compare_init(home,comp_dir)

# Plot:
ra.compare_mixed_traditional(home,comp_dir,tradpath_5,mixedpath_5,mymap,evaxlim,staxlim,paxlim,symbol_size,cbins)






#################################################################################
############## Comparison #1 - minimum of 30 stations per event #################
#################################################################################

# Model 1 - Inversion of Anza data
#   Residuals object:
tradpath_30=home+'v2anza2013_Mc8.5_pgrid_30sta_res4/v2anza2013_Mc8.5_pgrid_30sta_res4_robj.pckl'

# Model 2 - Mixed Effects inversion of Anza data
mixedpath_30=home+'mixedregr_v2anza2013_Mc_8.5_30sta/mixedregr_v2anza2013_30sta_Mc_8.5_VR_99.9_robj.pckl'

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
ra.compare_mixed_traditional(home,comp_dir,tradpath_30,mixedpath_30,mymap,evaxlim,staxlim,paxlim,symbol_size,cbins)




#################################################################################
############## Comparison #2 - minimum of 32 stations per event #################
#################################################################################
# Model 1 - Inversion of Anza data
#   Residuals object:
tradpath_32=home+'v2anza2013_Mc8.5_pgrid_32sta_res4/v2anza2013_Mc8.5_pgrid_32sta_res4_robj.pckl'

# Model 2 - Mixed Effects inversion of Anza data
mixedpath_32=home+'mixedregr_v2anza2013_Mc_8.5_32sta/mixedregr_v2anza2013_32sta_Mc_8.5_VR_99.9_robj.pckl'

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
ra.compare_mixed_traditional(home,comp_dir,tradpath_32,mixedpath_32,mymap,evaxlim,staxlim,paxlim,symbol_size,cbins)





##################################################################################
############# Comparison #3 - Between mixed models (5, 30, 32 sta) ###############
##################################################################################
##  How many models comparing?  will affect the plots...
#num_models_to_compare = 3
#
#tradpath_5=home+'v2anza2013_Mc8.5_pgrid_5sta_res4/v2anza2013_Mc8.5_pgrid_5sta_res4_robj.pckl'
#
## Model 2 - Mixed Effects inversion of Anza data
#mixedpath_5=home+'mixedregr_v2anza2013_Mc_8.5_res4/mixedregr_v2anza2013_Mc_8.5_VR_99.9_robj.pckl'
#
## Comparison directory name:
#comp_dir = 'compare_v2anza2013_5_30_32sta'
#
#################
###  Get info ##
#################
#
## Read in models:
#trfile=open(tradpath_5,'r')
#trobj_5sta=pickle.load(trfile)
#trfile.close()
#
#trfile=open(tradpath_30,'r')
#trobj_30sta=pickle.load(trfile)
#trfile.close()
#
#trfile=open(tradpath_32,'r')
#trobj_32sta=pickle.load(trfile)
#trfile.close()
#
######################
## Now mixed models:
#mifile=open(mixedpath_5,'r')
#miobj_5sta=pickle.load(mifile)
#mifile.close()
#
#mifile=open(mixedpath_30,'r')
#miobj_30sta=pickle.load(mifile)
#mifile.close()
#
#mifile=open(mixedpath_32,'r')
#miobj_32sta=pickle.load(mifile)
#mifile.close()
#
#####################
## Get number of random effects per db:
## 5 sta
#numev_persta_5 = ra.num_of_randeffects(trobj_5sta,numev=True)
#numsta_perev_5 = ra.num_of_randeffects(trobj_5sta,numsta=True)
#
## 30 sta
#numev_persta_30 = ra.num_of_randeffects(trobj_30sta,numev=True)
#numsta_perev_30 = ra.num_of_randeffects(trobj_30sta,numsta=True)
#
## 32 sta
#numev_persta_32 = ra.num_of_randeffects(trobj_32sta,numev=True)
#numsta_perev_32 = ra.num_of_randeffects(trobj_32sta,numsta=True)
#
#
#
##########################
######### Plotting #######
##########################
#
## First the events:
#
#
#mymap = 'Blues'
## 5sta, sinbgle-mean approach:
#nbins=[300,37]
#axlims=[[-4,4],[5,42]]
#clims=[-1000,2500]
#
#
#Z = [[0,0],[0,0]]
#levels = range(clims[0],clims[1]+1,1)
#CS3 = plt.contourf(Z, levels, cmap=mymap)
#plt.clf()
#
#evfig,evaxarr=plt.subplots(2,num_models_to_compare)
#
## Start plotting...start with single-mean models
#ra.compare_models_2dhist(np.unique(trobj_5sta.E_residual),numsta_perev_5,nbins,evaxarr[0,0],mymap,axlims,clims,'event','sta_per_ev')
#
#ra.compare_models_2dhist(np.unique(trobj_30sta.E_residual),numsta_perev_30,nbins,evaxarr[0,1],mymap,axlims,clims,'event','sta_per_ev')
#
#ra.compare_models_2dhist(np.unique(trobj_32sta.E_residual),numsta_perev_32,nbins,evaxarr[0,2],mymap,axlims,clims,'event','sta_per_ev')
#
## Now mixed effects:
#ra.compare_models_2dhist(np.unique(miobj_5sta.E_residual),numsta_perev_5,nbins,evaxarr[1,0],mymap,axlims,clims,'event','sta_per_ev')
#
#ra.compare_models_2dhist(np.unique(miobj_30sta.E_residual),numsta_perev_30,nbins,evaxarr[1,1],mymap,axlims,clims,'event','sta_per_ev')
#
#ra.compare_models_2dhist(np.unique(miobj_32sta.E_residual),numsta_perev_32,nbins,evaxarr[1,2],mymap,axlims,clims,'event','sta_per_ev')
#
#
#evfig.colorbar(CS3)
#
#evfig.title('Event Terms - Single-mean Model (top), Mixed Effects Model (bottom)')