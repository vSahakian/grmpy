#
#
# Compare event, site, and path terms between models
#   1.  Plot 2d histograms of event and site term, for the residual vs. num ev/sta or sta/ev, compare db's
#           - In one plot, top row: single-mean.  bottom row: mixed.  columns: 5 sta, 30 sta, 32 sta
#   2.  Plot 1d pdf's and histograms of event, path and site terms between db's for single-mean models (left), and mixed (right).
#           - One plot per term (event, path, site)
#   3.  Plot 1d pdf's and histograms for just the 5sta db
#           - Three plots, one per term, plot single-mean and mixed on the same plot for each db
#   4.  Plot just the residuals between the same models, for different db's 
#             (i.e., single (5) vs. single (30))
#           - Makes three plots for each time two models are compared
#
#






import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import res_analysis as ra
from scipy.stats import norm
#from statsmodels.stats.weightstats.CompareMeans import ztest_ind

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

# Model 1 - Inversion of Anza data with 30 statoins per event
#   Single-mean Residuals object:
tradpath_30=home+'v2anza2013_Mc8.5_pgrid_30sta_res4/v2anza2013_Mc8.5_pgrid_30sta_res4_robj.pckl'

# Model 2 - Mixed Effects inversion of Anza data
mixedpath_30=home+'mixedregr_v2anza2013_Mc_8.5_30sta/mixedregr_v2anza2013_30sta_Mc_8.5_VR_99.9_robj.pckl'


# Model 1 - Inversion of Anza data with 32 stations per event
#   Single-mean Residuals object:
tradpath_32=home+'v2anza2013_Mc8.5_pgrid_32sta_res4/v2anza2013_Mc8.5_pgrid_32sta_res4_robj.pckl'

# Model 2 - Mixed Effects inversion of Anza data
mixedpath_32=home+'mixedregr_v2anza2013_Mc_8.5_32sta/mixedregr_v2anza2013_32sta_Mc_8.5_VR_99.9_robj.pckl'


#################################################################################
############ Comparison #3 - Between mixed models (5, 30, 32 sta) ###############
#################################################################################

# Model with 5 stations per event
# Single-mean residuals object
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
#    Such as:
#        Number of events per station (for each database), 
#        Number of stations recording each event (for each database)
#        Unique events/stations for these
# 
# For database:
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


###################################
######## Initiate directory #######
###################################

# Initiate directory:
ra.compare_init(home,comp_dir)



#########################
######## Plotting #######
#########################

print 'plotting event histograms'
# Initiate directory:
ra.compare_init(home,comp_dir)


#### First the events ####

# Plotting parameters:
mymap = 'gnuplot'
# 5sta, sinbgle-mean approach:
# nuber of bins, in x and y
nbins=[300,37]
# Axis limits [[resmin, resmax],[num_stations_min,num_stations_ymax]]
axlims=[[-4,4],[5,42]]
# Color limits (for histogram)
clims=[-1000,2500]

# figure size:
fsize=(16,13)


## Make "fake plot" for colorbar
Z = [[0,0],[0,0]]
levels = range(clims[0],clims[1]+1,1)
CS3 = plt.contourf(Z, levels, cmap=mymap)
plt.clf()

# Make the subplots for the events/initiate the figure
evfig,evaxarr=plt.subplots(2,num_models_to_compare,figsize=fsize)

# Start plotting...start with single-mean models
# plot event residuals/number of stations per vent for the 5 sta database
ra.compare_models_2dhist(np.array(trobj_5sta.E_residual)[ev5unique],numsta_perev_5,nbins,evaxarr[0,0],mymap,axlims,clims,'event','sta_per_ev')
# for the 30 sta db
ra.compare_models_2dhist(np.array(trobj_30sta.E_residual)[ev30unique],numsta_perev_30,nbins,evaxarr[0,1],mymap,axlims,clims,'event','sta_per_ev')
# for the 32 sta db
ra.compare_models_2dhist(np.array(trobj_32sta.E_residual)[ev32unique],numsta_perev_32,nbins,evaxarr[0,2],mymap,axlims,clims,'event','sta_per_ev')

# Now mixed effects (same as above):
ra.compare_models_2dhist(np.array(miobj_5sta.E_residual)[ev5unique],numsta_perev_5,nbins,evaxarr[1,0],mymap,axlims,clims,'event','sta_per_ev')

ra.compare_models_2dhist(np.array(miobj_30sta.E_residual)[ev30unique],numsta_perev_30,nbins,evaxarr[1,1],mymap,axlims,clims,'event','sta_per_ev')

ra.compare_models_2dhist(np.array(miobj_32sta.E_residual)[ev32unique],numsta_perev_32,nbins,evaxarr[1,2],mymap,axlims,clims,'event','sta_per_ev')

# Add the colorbar
evfig.colorbar(CS3)

#evfig.title('Event Terms - Single-mean Model (top), Mixed Effects Model (bottom)')
# Save the figure (png and pdf)
evfig.savefig(home+comp_dir+'/figs/events.png')
evfig.savefig(home+comp_dir+'/figs/pdfs/events.pdf')

#############################################################3
print 'plotting site histograms'
# Then site histograms:

# Plotting parameters 
mymap = 'gnu_plot'
# 5sta, sinbgle-mean approach:
nbins=[8,8190]
axlims=[[-4,4],[2,5000]]
clims=[0,2]
# Figure size
fsize=(16,13)

# Make fake plot for the colorbar
Z = [[0,0],[0,0]]
levels = range(clims[0],clims[1]+1,1)
CS3 = plt.contourf(Z, levels, cmap=mymap)
plt.clf()

# Initiate the figure/subplots
stfig,staxarr=plt.subplots(2,num_models_to_compare,figsize=fsize)

# Start plotting...start with single-mean models
ra.compare_models_2dhist(np.array(trobj_5sta.site_terms)[st5unique],numev_persta_5,nbins,staxarr[0,0],mymap,axlims,clims,'site','ev_per_sta')

ra.compare_models_2dhist(np.array(trobj_30sta.site_terms)[st30unique],numev_persta_30,nbins,staxarr[0,1],mymap,axlims,clims,'site','ev_per_sta')

ra.compare_models_2dhist(np.array(trobj_32sta.site_terms)[st32unique],numev_persta_32,nbins,staxarr[0,2],mymap,axlims,clims,'site','ev_per_sta')

# Now mixed effects:
ra.compare_models_2dhist(np.array(miobj_5sta.site_terms)[st5unique],numev_persta_5,nbins,staxarr[1,0],mymap,axlims,clims,'site','ev_per_sta')

ra.compare_models_2dhist(np.array(miobj_30sta.site_terms)[st30unique],numev_persta_30,nbins,staxarr[1,1],mymap,axlims,clims,'site','ev_per_sta')

ra.compare_models_2dhist(np.array(miobj_32sta.site_terms)[st32unique],numev_persta_32,nbins,staxarr[1,2],mymap,axlims,clims,'site','ev_per_sta')

stfig.colorbar(CS3)


#evfig.title('Event Terms - Single-mean Model (top), Mixed Effects Model (bottom)')
# Save figure to png and pdf
stfig.savefig(home+comp_dir+'/figs/sites.png')
stfig.savefig(home+comp_dir+'/figs/pdfs/sites.pdf')



#############################################################3
##########      Make plots of 1-d pdf/histogram   ###########3
#############################################################3



# First path....
#   Can't plot 2d histograms for path since don't have "number of stations per path", etc..
print 'plotting pdfs and 1d histograms of path'

# Getting paramters for 5 sta per event db
# First get sorted values of path term, first for the single-mean(tr):
p_tr_sorted_5sta = np.sort(trobj_5sta.path_terms)
# Then for mixed
p_mi_sorted_5sta = np.sort(miobj_5sta.path_terms)
# Also get the pdf of the single-ean and mixed
ptr5pdf = norm.pdf(p_tr_sorted_5sta)
pmi5pdf = norm.pdf(p_mi_sorted_5sta)

# Same paramters for 30 sta/event db
p_tr_sorted_30sta = np.sort(trobj_30sta.path_terms)
p_mi_sorted_30sta = np.sort(miobj_30sta.path_terms)
ptr30pdf = norm.pdf(p_tr_sorted_30sta)
pmi30pdf = norm.pdf(p_mi_sorted_30sta)

# Same parameters for 32 sta/event db
p_tr_sorted_32sta = np.sort(trobj_32sta.path_terms)
p_mi_sorted_32sta = np.sort(miobj_32sta.path_terms)
ptr32pdf = norm.pdf(p_tr_sorted_32sta)
pmi32pdf = norm.pdf(p_mi_sorted_32sta)

#####################3
#### Plotting ########
######################

# Initiate figure; left for single-mean, right for mixed
fsize=(20,13)
# get figure with subplots...
pfig,axarr = plt.subplots(1,2,figsize=fsize)

### SINGLE-MEAN, LEFT COLUMN ###
# First plot the pdf and histogram of the 5sta, on the left subplot:
tr5sta=axarr[0].plot(p_tr_sorted_5sta,ptr5pdf,linestyle='-',color='gray',label='5sta')
tr30sta=axarr[0].plot(p_tr_sorted_30sta,ptr30pdf,linestyle='--',color='red',label='30sta')
tr32sta=axarr[0].plot(p_tr_sorted_32sta,ptr32pdf,linestyle='-.',color='blue',label='32sta')

# also on this left subplot, plot the histogram beneath the pdf
axarr[0].hist(p_tr_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='gray',normed=True)
axarr[0].hist(p_tr_sorted_30sta,bins=50,alpha=0.2,histtype='stepfilled',color='red',normed=True)
axarr[0].hist(p_tr_sorted_32sta,bins=50,alpha=0.2,histtype='stepfilled',color='blue',normed=True)

# Get legend labels
axhandle,axlabel=axarr[0].get_legend_handles_labels()
# Make legend
axarr[0].legend(axhandle,axlabel,loc='upper right',frameon=False)

# Label plot
axarr[0].set_xlabel('Path Term (ln residual')
axarr[0].set_ylabel('Probability Density/Normed count')
axarr[0].title.set_text('Path PDF and Normed histogram for Single-mean model')

### MIXED MODEL, RIGHT COLUMN ###
# First plot the pdf of the 5sta:
mi5sta=axarr[1].plot(p_mi_sorted_5sta,pmi5pdf,linestyle='-',color='gray',label='5sta')
# Then the 30 sta
mi30sta=axarr[1].plot(p_mi_sorted_30sta,pmi30pdf,linestyle='--',color='red',label='30sta')
# Then the 32 sta 
mi32sta=axarr[1].plot(p_mi_sorted_32sta,pmi32pdf,linestyle='-.',color='blue',label='32sta')

# Now plot histograms for 5, 30, and 32 sta db's on same subplot 
axarr[1].hist(p_mi_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='gray',normed=True)
axarr[1].hist(p_mi_sorted_30sta,bins=50,alpha=0.2,histtype='stepfilled',color='red',normed=True)
axarr[1].hist(p_mi_sorted_32sta,bins=50,alpha=0.2,histtype='stepfilled',color='blue',normed=True)

# Get axis legend data...
axhandle,axlabel=axarr[0].get_legend_handles_labels()
# Make legend
axarr[1].legend(axhandle,axlabel,loc='upper right',frameon=False)

# Set labels
axarr[1].set_xlabel('Path Term (ln residual')
axarr[1].set_ylabel('Probability Density/Normed count')
axarr[1].title.set_text('Path PDF and Normed histogram for Mixed effects model')

### Save figure ###
pfig.savefig(home+comp_dir+'/figs/path_pdf.png')
pfig.savefig(home+comp_dir+'/figs/pdfs/path_pdf.pdf')


#############################################################3
#############################################################3
#############################################################3
# Make Same figure for stations:
print 'plotting pdfs and 1d histograms of site'

# 5 sta db
# Sort the single-mean and mixed model site terms
s_tr_sorted_5sta = np.sort(trobj_5sta.site_terms)
s_mi_sorted_5sta = np.sort(miobj_5sta.site_terms)
# Get the pdf's of these 
str5pdf = norm.pdf(s_tr_sorted_5sta)
smi5pdf = norm.pdf(s_mi_sorted_5sta)

# 30 sta db
# Sort the single-mean and mixed model site terms
s_tr_sorted_30sta = np.sort(trobj_30sta.site_terms)
s_mi_sorted_30sta = np.sort(miobj_30sta.site_terms)
# Get the pdf's of these 
str30pdf = norm.pdf(s_tr_sorted_30sta)
smi30pdf = norm.pdf(s_mi_sorted_30sta)

# 32 sta db
# Sort the single-mean and mixed model site terms
s_tr_sorted_32sta = np.sort(trobj_32sta.site_terms)
s_mi_sorted_32sta = np.sort(miobj_32sta.site_terms)
# Get the pdf's of these 
str32pdf = norm.pdf(s_tr_sorted_32sta)
smi32pdf = norm.pdf(s_mi_sorted_32sta)

#####################3
#### Plotting ########
######################

# Initiate figure; left for single-mean, right for mixed
fsize=(20,13)
sfig,axarr = plt.subplots(1,2,figsize=fsize)

# First plot the pdf of the single-mean models:
# 5sta
tr5sta=axarr[0].plot(s_tr_sorted_5sta,str5pdf,linestyle='-',color='gray',label='5sta')
# 30 sta
tr30sta=axarr[0].plot(s_tr_sorted_30sta,str30pdf,linestyle='--',color='red',label='30sta')
# 32 sta
tr32sta=axarr[0].plot(s_tr_sorted_32sta,str32pdf,linestyle='-.',color='blue',label='32sta')

# Then on teh samme plot, plot the histograms for each database (5sta, 30sta, 32sta)
axarr[0].hist(s_tr_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='gray',normed=True)
axarr[0].hist(s_tr_sorted_30sta,bins=50,alpha=0.2,histtype='stepfilled',color='red',normed=True)
axarr[0].hist(s_tr_sorted_32sta,bins=50,alpha=0.2,histtype='stepfilled',color='blue',normed=True)

# Get legend labels
axhandle,axlabel=axarr[0].get_legend_handles_labels()
# Make legend
axarr[0].legend(axhandle,axlabel,loc='upper right',frameon=False)

# Label subplot
axarr[0].set_xlabel('Site Term (ln residual')
axarr[0].set_ylabel('Probability Density/Normed count')
axarr[0].title.set_text('Site PDF and Normed histogram for Single-mean model')

#############
#  Now for mixed:
# First plot the pdf of the mixed model 5sta, 30 sta, 32 sta:
mi5sta=axarr[1].plot(s_mi_sorted_5sta,smi5pdf,linestyle='-',color='gray',label='5sta')
mi30sta=axarr[1].plot(s_mi_sorted_30sta,smi30pdf,linestyle='--',color='red',label='30sta')
mi32sta=axarr[1].plot(s_mi_sorted_32sta,smi32pdf,linestyle='-.',color='blue',label='32sta')

# then the histograms of the mixed model 5 sta, 30 sta, and 32 sta
axarr[1].hist(s_mi_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='gray',normed=True)
axarr[1].hist(s_mi_sorted_30sta,bins=50,alpha=0.2,histtype='stepfilled',color='red',normed=True)
axarr[1].hist(s_mi_sorted_32sta,bins=50,alpha=0.2,histtype='stepfilled',color='blue',normed=True)

# Get axis legend info
axhandle,axlabel=axarr[0].get_legend_handles_labels()
# Plot legend
axarr[1].legend(axhandle,axlabel,loc='upper right',frameon=False)

# Plot labels
axarr[1].set_xlabel('Site Term (ln residual')
axarr[1].set_ylabel('Probability Density/Normed count')
axarr[1].title.set_text('Site PDF and Normed histogram for Mixed effects model')

### Save ###
sfig.savefig(home+comp_dir+'/figs/site_pdf.png')
sfig.savefig(home+comp_dir+'/figs/pdfs/site_pdf.png')

#############################################################3
#############################################################3
#############################################################3
# Make Same figure for events:
print 'plotting pdfs and 1d histograms of events'

# Get sorted values of event term for slingle-mean and mixed
# Also get pdf...

# sorted 5sta...
e_tr_sorted_5sta = np.sort(trobj_5sta.E_residual)
e_mi_sorted_5sta = np.sort(miobj_5sta.E_residual)
# pdf 5 sta
etr5pdf = norm.pdf(e_tr_sorted_5sta)
emi5pdf = norm.pdf(e_mi_sorted_5sta)

# Sorted 30 sta
e_tr_sorted_30sta = np.sort(trobj_30sta.E_residual)
e_mi_sorted_30sta = np.sort(miobj_30sta.E_residual)
# pdf 30 sta
etr30pdf = norm.pdf(e_tr_sorted_30sta)
emi30pdf = norm.pdf(e_mi_sorted_30sta)

# Sorted 32 sta
e_tr_sorted_32sta = np.sort(trobj_32sta.E_residual)
e_mi_sorted_32sta = np.sort(miobj_32sta.E_residual)
# pdf 32 sta
etr32pdf = norm.pdf(e_tr_sorted_32sta)
emi32pdf = norm.pdf(e_mi_sorted_32sta)

#####################3
#### Plotting ########
######################

# Initiate figure; left for single-mean, right for mixed
fsize=(20,13)
efig,axarr = plt.subplots(1,2,figsize=fsize)

#####
# Single mean plots, left subplot:

# First plot the pdf of the 5sta, 30 sta, 32 sta:
tr5sta=axarr[0].plot(e_tr_sorted_5sta,etr5pdf,linestyle='-',color='gray',label='5sta')
tr30sta=axarr[0].plot(e_tr_sorted_30sta,etr30pdf,linestyle='--',color='red',label='30sta')
tr32sta=axarr[0].plot(e_tr_sorted_32sta,etr32pdf,linestyle='-.',color='blue',label='32sta')

# Then the 1d histogram of the single-mean 5 sta, 30 sta, and 32 sta
axarr[0].hist(e_tr_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='gray',normed=True)
axarr[0].hist(e_tr_sorted_30sta,bins=50,alpha=0.2,histtype='stepfilled',color='red',normed=True)
axarr[0].hist(e_tr_sorted_32sta,bins=50,alpha=0.2,histtype='stepfilled',color='blue',normed=True)

# Get axis legend data, plto legend:
axhandle,axlabel=axarr[0].get_legend_handles_labels()
axarr[0].legend(axhandle,axlabel,loc='upper right',frameon=False)

# Set axis labels
axarr[0].set_xlabel('Event Term (ln residual')
axarr[0].set_ylabel('Probability Density/Normed count')
axarr[0].title.set_text('Event PDF and Normed histogram for Single-mean model')

#############
#  Mixed plots, right subplot:

# First plot the pdf of the mixed models for 5sta, 30 sta, 32 sta:
mi5sta=axarr[1].plot(e_mi_sorted_5sta,emi5pdf,linestyle='-',color='gray',label='5sta')
mi30sta=axarr[1].plot(e_mi_sorted_30sta,emi30pdf,linestyle='--',color='red',label='30sta')
mi32sta=axarr[1].plot(e_mi_sorted_32sta,emi32pdf,linestyle='-.',color='blue',label='32sta')

# Then the 1d histograms for the mixed models, 5sta, 30 sta, 32 sta db's
axarr[1].hist(e_mi_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='gray',normed=True)
axarr[1].hist(e_mi_sorted_30sta,bins=50,alpha=0.2,histtype='stepfilled',color='red',normed=True)
axarr[1].hist(e_mi_sorted_32sta,bins=50,alpha=0.2,histtype='stepfilled',color='blue',normed=True)

# Get legend handles, plot legend
axhandle,axlabel=axarr[0].get_legend_handles_labels()
axarr[1].legend(axhandle,axlabel,loc='upper right',frameon=False)

# Plot labels
axarr[1].set_xlabel('Event Term (ln residual')
axarr[1].set_ylabel('Probability Density/Normed count')
axarr[1].title.set_text('Event PDF and Normed histogram for Mixed effects model')

### Save ###
efig.savefig(home+comp_dir+'/figs/event_pdf.png')
efig.savefig(home+comp_dir+'/figs/pdfs/event_pdf.pdf')



#############################################################3
############## 5 sta Single mean and Mixed pdfs ##############
#############################################################3

# Make a single figure (one for event, one for site, one for path) to compare mixed vs. single-mean
#    between just the 5 sta case:

# Events:

# Initiate figure
ev5stafig=plt.figure()

# plot pdfs
etr5sta=plt.plot(e_tr_sorted_5sta,etr5pdf,linestyle='--',color='red',label='single-mean')
emi5sta=plt.plot(e_mi_sorted_5sta,emi5pdf,linestyle='-.',color='blue',label='mixed')
# plot hists
plt.hist(e_tr_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='red',normed=True)
plt.hist(e_mi_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='blue',normed=True)

# Get z test value???
e5sta_z = ra.z_test(e_tr_sorted_5sta,e_mi_sorted_5sta)

# Plot legend
plt.legend(loc='upper right',frameon=False)

# Plot axis labels
plt.xlabel('Event Term (ln residual)')
plt.ylabel('Probability/Normed count')
plt.title('Event Terms: Single-mean vs. Mixed effects') # \n Z-test: %f'  % e5sta_z)

# Save
ev5stafig.savefig(home+comp_dir+'/figs/event_5sta.png')
ev5stafig.savefig(home+comp_dir+'/figs/pdfs/event_5sta.pdf')

######
# Stations:

# Initiate plot
st5stafig=plt.figure()

# plot pdfs (single-mean, mixed)
str5sta=plt.plot(s_tr_sorted_5sta,str5pdf,linestyle='--',color='red',label='single-mean')
smi5sta=plt.plot(s_mi_sorted_5sta,smi5pdf,linestyle='-.',color='blue',label='mixed')
# plot hists
plt.hist(s_tr_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='red',normed=True)
plt.hist(s_mi_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='blue',normed=True)

# Get z test value???
s5sta_z  = ra.z_test(s_tr_sorted_5sta,s_mi_sorted_5sta)

# Plot legend
plt.legend(loc='upper right',frameon=False)

# plot axis labels
plt.xlabel('Site Term (ln residual)')
plt.ylabel('Probability/Normed count')
plt.title('Site Terms: Single-mean vs. Mixed effects') # \n Z-test: %f' % s5sta_z)

# Save site figure
st5stafig.savefig(home+comp_dir+'/figs/site_5sta.png')
st5stafig.savefig(home+comp_dir+'/figs/pdfs/site_5sta.pdf')

######
# Paths:

# Initiate figure
pa5stafig=plt.figure()

# plot pdfs
ptr5sta=plt.plot(p_tr_sorted_5sta,ptr5pdf,linestyle='--',color='red',label='single-mean')
pmi5sta=plt.plot(p_mi_sorted_5sta,pmi5pdf,linestyle='-.',color='blue',label='mixed')
# plot hists
plt.hist(p_tr_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='red',normed=True)
plt.hist(p_mi_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='blue',normed=True)

# Get z test value???
p5sta_z = ra.z_test(p_tr_sorted_5sta,p_mi_sorted_5sta)

# Plot legend
plt.legend(loc='upper right',frameon=False)

# Plot labels
plt.xlabel('Path Term (ln residual)')
plt.ylabel('Probability/Normed count')
plt.title('Path Terms: Single-mean vs. Mixed effects') # \n Z-test: %f' % p5sta_z)

# Save
pa5stafig.savefig(home+comp_dir+'/figs/path_5sta.png')
pa5stafig.savefig(home+comp_dir+'/figs/pdfs/path_5sta.pdf')



#################################################################################
#################################################################################
############# Comparison - Single vs. single model for 5, 30, 32 ################
#################################################################################
#################################################################################


#   1.  Compare single-mean vs. single-mean for 5 and 30 sta
method_type='single-mean'
mymap='gnuplot'
# Axis limits for event terms, station temrs, and path terms:
evaxlim=[[-2,2],[-2,2]]
staaxlim=[[-2.5,2.5],[-2.5,2.5]]
pathaxlim=[[-3,3],[-3,3]]
# Symbol size for event, path, site
symbol_size=[20,10,20]
# Bin size for colorscale for event, site
cbins=[.01,.1]
# limits for colorscale [[cmin_numsta_per_event,cmax_numsta_per_event],[cmin_numev_per_station,cmax_numev_per_station]]:
clims=[[30,37],[1,70]]


# Get the indices where the smaller databse is in the larger databse to plot together
comp_indices = ra.get_comparison_indices(trobj_5sta,trobj_30sta)
# Plot single vs. single first....
ra.plot_method_method(home,comp_dir,method_type,tradpath_5,tradpath_30,'5sta','30sta',comp_indices,mymap,evaxlim,staaxlim,pathaxlim,symbol_size,cbins,clims)

# Now plot mixed:
method_type='mixed'
ra.plot_method_method(home,comp_dir,method_type,mixedpath_5,mixedpath_30,'5sta','30sta',comp_indices,mymap,evaxlim,staaxlim,pathaxlim,symbol_size,cbins,clims)


##################################
##################################
##################################


#   2.  Compare single-mean vs. single-mean for 5 and 32 sta
method_type='single-mean'
mymap='gnuplot'
# Axis limits for event terms, station temrs, and path terms:
evaxlim=[[-2,2],[-2,2]]
staaxlim=[[-2.5,2.5],[-2.5,2.5]]
pathaxlim=[[-3,3],[-3,3]]
# Symbol size for event, path, site
symbol_size=[20,10,20]
# Bin size for colorscale for event, site
cbins=[.01,.011]
# limits for colorscale [[cmin_numsta_per_event,cmax_numsta_per_event],[cmin_numev_per_station,cmax_numev_per_station]]:
clims=[[32,37],[1,25]]


# Get the indices where the smaller databse is in the larger databse to plot together
comp_indices = ra.get_comparison_indices(trobj_5sta,trobj_32sta)
# Plot single vs. single first....
ra.plot_method_method(home,comp_dir,method_type,tradpath_5,tradpath_32,'5sta','32sta',comp_indices,mymap,evaxlim,staaxlim,pathaxlim,symbol_size,cbins,clims)

# Now plot mixed:
method_type='mixed'
ra.plot_method_method(home,comp_dir,method_type,mixedpath_5,mixedpath_32,'5sta','32sta',comp_indices,mymap,evaxlim,staaxlim,pathaxlim,symbol_size,cbins,clims)



##################################
##################################
##################################


#   3.  Compare single-mean vs. single-mean for 30 and 32 sta
method_type='single-mean'
mymap='gnuplot'
# Axis limits for event terms, station temrs, and path terms:
evaxlim=[[-2,2],[-2,2]]
staaxlim=[[-2.5,2.5],[-2.5,2.5]]
pathaxlim=[[-3,3],[-3,3]]
# Symbol size for event, path, site
symbol_size=[20,10,20]
# Bin size for colorscale for event, site
cbins=[.01,.1]
# limits for colorscale [[cmin_numsta_per_event,cmax_numsta_per_event],[cmin_numev_per_station,cmax_numev_per_station]]:
clims=[[32,35],[1,20]]


# Get the indices where the smaller databse is in the larger databse to plot together
comp_indices = ra.get_comparison_indices(trobj_30sta,trobj_32sta)
# Plot single vs. single first....
ra.plot_method_method(home,comp_dir,method_type,tradpath_30,tradpath_32,'30sta','32sta',comp_indices,mymap,evaxlim,staaxlim,pathaxlim,symbol_size,cbins,clims)

# Now plot mixed:
method_type='mixed'
ra.plot_method_method(home,comp_dir,method_type,mixedpath_30,mixedpath_32,'30sta','32sta',comp_indices,mymap,evaxlim,staaxlim,pathaxlim,symbol_size,cbins,clims)
