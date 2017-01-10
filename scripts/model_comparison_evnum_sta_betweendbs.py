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

#################################################################################
############## Comparison #1 - minimum of 30 stations per event #################
#################################################################################

# Model 1 - Inversion of Anza data
#   Residuals object:
tradpath_30=home+'v2anza2013_Mc8.5_pgrid_30sta_res4/v2anza2013_Mc8.5_pgrid_30sta_res4_robj.pckl'

# Model 2 - Mixed Effects inversion of Anza data
mixedpath_30=home+'mixedregr_v2anza2013_Mc_8.5_30sta/mixedregr_v2anza2013_30sta_Mc_8.5_VR_99.9_robj.pckl'


# Model 1 - Inversion of Anza data
#   Residuals object:
tradpath_32=home+'v2anza2013_Mc8.5_pgrid_32sta_res4/v2anza2013_Mc8.5_pgrid_32sta_res4_robj.pckl'

# Model 2 - Mixed Effects inversion of Anza data
mixedpath_32=home+'mixedregr_v2anza2013_Mc_8.5_32sta/mixedregr_v2anza2013_32sta_Mc_8.5_VR_99.9_robj.pckl'


#################################################################################
############ Comparison #3 - Between mixed models (5, 30, 32 sta) ###############
#################################################################################

tradpath_5=home+'v2anza2013_Mc8.5_pgrid_5sta_res4/v2anza2013_Mc8.5_pgrid_5sta_res4_robj.pckl'

# Model 2 - Mixed Effects inversion of Anza data
mixedpath_5=home+'mixedregr_v2anza2013_Mc_8.5_res4/mixedregr_v2anza2013_Mc_8.5_VR_99.9_robj.pckl'


#  How many models comparing?  will affect the plots...
num_models_to_compare = 3
# Comparison directory name:
comp_dir='compare_v2anza2013_5_30_32'


################
##  Get info ##
################

print 'Reading in models...'
# Read in models:
trfile=open(tradpath_5,'r')
trobj_5sta=pickle.load(trfile)
trfile.close()

trfile=open(tradpath_30,'r')
trobj_30sta=pickle.load(trfile)
trfile.close()

trfile=open(tradpath_32,'r')
trobj_32sta=pickle.load(trfile)
trfile.close()

#####################
print 'Reading in mixed models...'
# Now mixed models:
mifile=open(mixedpath_5,'r')
miobj_5sta=pickle.load(mifile)
mifile.close()

mifile=open(mixedpath_30,'r')
miobj_30sta=pickle.load(mifile)
mifile.close()

mifile=open(mixedpath_32,'r')
miobj_32sta=pickle.load(mifile)
mifile.close()

####################
print 'getting aux info'
# Get number of random effects per db:
# 5 sta
numev_persta_5 = ra.num_of_randeffects(trobj_5sta,numev=True)
numsta_perev_5 = ra.num_of_randeffects(trobj_5sta,numsta=True)
unev,ev5unique = np.unique(trobj_5sta.evnum,return_index=True)
unst,st5unique = np.unique(trobj_5sta.sta,return_index=True)

# 30 sta
numev_persta_30 = ra.num_of_randeffects(trobj_30sta,numev=True)
numsta_perev_30 = ra.num_of_randeffects(trobj_30sta,numsta=True)
unev,ev30unique = np.unique(trobj_30sta.evnum,return_index=True)
unst,st30unique = np.unique(trobj_30sta.sta,return_index=True)

# 32 sta
numev_persta_32 = ra.num_of_randeffects(trobj_32sta,numev=True)
numsta_perev_32 = ra.num_of_randeffects(trobj_32sta,numsta=True)
unev,ev32unique = np.unique(trobj_32sta.evnum,return_index=True)
unst,st32unique = np.unique(trobj_32sta.sta,return_index=True)


#########################
######## Plotting #######
#########################

print 'plotting events'
# Initiate directory:
ra.compare_init(home,comp_dir)


# First the events:
mymap = 'Blues'
# 5sta, sinbgle-mean approach:
nbins=[300,37]
axlims=[[-4,4],[5,42]]
clims=[-1000,2500]

fsize=(16,13)

Z = [[0,0],[0,0]]
levels = range(clims[0],clims[1]+1,1)
CS3 = plt.contourf(Z, levels, cmap=mymap)
plt.clf()

evfig,evaxarr=plt.subplots(2,num_models_to_compare,figsize=fsize)

# Start plotting...start with single-mean models
ra.compare_models_2dhist(np.array(trobj_5sta.E_residual)[ev5unique],numsta_perev_5,nbins,evaxarr[0,0],mymap,axlims,clims,'event','sta_per_ev')

ra.compare_models_2dhist(np.array(trobj_30sta.E_residual)[ev30unique],numsta_perev_30,nbins,evaxarr[0,1],mymap,axlims,clims,'event','sta_per_ev')

ra.compare_models_2dhist(np.array(trobj_32sta.E_residual)[ev32unique],numsta_perev_32,nbins,evaxarr[0,2],mymap,axlims,clims,'event','sta_per_ev')

# Now mixed effects:
ra.compare_models_2dhist(np.array(miobj_5sta.E_residual)[ev5unique],numsta_perev_5,nbins,evaxarr[1,0],mymap,axlims,clims,'event','sta_per_ev')

ra.compare_models_2dhist(np.array(miobj_30sta.E_residual)[ev30unique],numsta_perev_30,nbins,evaxarr[1,1],mymap,axlims,clims,'event','sta_per_ev')

ra.compare_models_2dhist(np.array(miobj_32sta.E_residual)[ev32unique],numsta_perev_32,nbins,evaxarr[1,2],mymap,axlims,clims,'event','sta_per_ev')


evfig.colorbar(CS3)

#evfig.title('Event Terms - Single-mean Model (top), Mixed Effects Model (bottom)')
evfig.savefig(home+comp_dir+'/figs/events.png')

#############################################################3
print 'plotting sites'
# Then sites:
mymap = 'Blues'
# 5sta, sinbgle-mean approach:
nbins=[300,8190]
axlims=[[-4,4],[2,8192]]
clims=[-1000,2500]

fsize=(16,13)

Z = [[0,0],[0,0]]
levels = range(clims[0],clims[1]+1,1)
CS3 = plt.contourf(Z, levels, cmap=mymap)
plt.clf()

stfig,staxarr=plt.subplots(2,num_models_to_compare,figsize=fsize)

# Start plotting...start with single-mean models
ra.compare_models_2dhist(np.array(trobj_5sta.site_term)[ev5unique],numev_persta_5,nbins,staxarr[0,0],mymap,axlims,clims,'site','ev_per_sta')

ra.compare_models_2dhist(np.array(trobj_30sta.site_term)[ev30unique],numev_persta_30,nbins,staxarr[0,1],mymap,axlims,clims,'site','ev_per_sta')

ra.compare_models_2dhist(np.array(trobj_32sta.site_term)[ev32unique],numev_persta_32,nbins,staxarr[0,2],mymap,axlims,clims,'site','ev_per_sta')

# Now mixed effects:
ra.compare_models_2dhist(np.array(miobj_5sta.site_term)[ev5unique],numev_persta_5,nbins,staxarr[1,0],mymap,axlims,clims,'site','ev_per_sta')

ra.compare_models_2dhist(np.array(miobj_30sta.site_term)[ev30unique],numev_persta_30,nbins,staxarr[1,1],mymap,axlims,clims,'site','ev_per_sta')

ra.compare_models_2dhist(np.array(miobj_32sta.site_term)[ev32unique],numev_persta_32,nbins,staxarr[1,2],mymap,axlims,clims,'site','ev_per_sta')

#evfig.title('Event Terms - Single-mean Model (top), Mixed Effects Model (bottom)')
stfig.savefig(home+comp_dir+'/figs/sites.png')