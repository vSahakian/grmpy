# Compare event, site, and path terms between models
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
nbins=[8,8190]
axlims=[[-4,4],[2,5000]]
clims=[0,2]

fsize=(16,13)

Z = [[0,0],[0,0]]
levels = range(clims[0],clims[1]+1,1)
CS3 = plt.contourf(Z, levels, cmap=mymap)
plt.clf()

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
stfig.savefig(home+comp_dir+'/figs/sites.png')


#############################################################3
#############################################################3
#############################################################3
# Now path....

# First get sorted values of path term:
p_tr_sorted_5sta = np.sort(trobj_5sta.path_terms)
p_mi_sorted_5sta = np.sort(miobj_5sta.path_terms)
ptr5pdf = norm.pdf(p_tr_sorted_5sta)
pmi5pdf = norm.pdf(p_mi_sorted_5sta)

p_tr_sorted_30sta = np.sort(trobj_30sta.path_terms)
p_mi_sorted_30sta = np.sort(miobj_30sta.path_terms)
ptr30pdf = norm.pdf(p_tr_sorted_30sta)
pmi30pdf = norm.pdf(p_mi_sorted_30sta)

p_tr_sorted_32sta = np.sort(trobj_32sta.path_terms)
p_mi_sorted_32sta = np.sort(miobj_32sta.path_terms)
ptr32pdf = norm.pdf(p_tr_sorted_32sta)
pmi32pdf = norm.pdf(p_mi_sorted_32sta)

#####################3
#### Plotting ########
######################

# Initiate figure; left for single-mean, right for mixed
fsize=(20,13)
pfig,axarr = plt.subplots(1,2,figsize=fsize)

# First plot the pdf and histogram of the 5sta:
tr5sta=axarr[0].plot(p_tr_sorted_5sta,ptr5pdf,linestyle='-',color='gray',label='5sta')
tr30sta=axarr[0].plot(p_tr_sorted_30sta,ptr30pdf,linestyle='--',color='red',label='30sta')
tr32sta=axarr[0].plot(p_tr_sorted_32sta,ptr32pdf,linestyle='-.',color='blue',label='32sta')

axarr[0].hist(p_tr_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='gray',normed=True)
axarr[0].hist(p_tr_sorted_30sta,bins=50,alpha=0.2,histtype='stepfilled',color='red',normed=True)
axarr[0].hist(p_tr_sorted_32sta,bins=50,alpha=0.2,histtype='stepfilled',color='blue',normed=True)

axhandle,axlabel=axarr[0].get_legend_handles_labels()
axarr[0].legend(axhandle,axlabel,loc='upper right',frameon=False)

axarr[0].set_xlabel('Path Term (ln residual')
axarr[0].set_ylabel('Probability Density/Normed count')
axarr[0].title.set_text('Path PDF and Normed histogram for Single-mean model')

#############
#  Now for mixed:
# First plot the pdf and histogram of the 5sta:
mi5sta=axarr[1].plot(p_mi_sorted_5sta,pmi5pdf,linestyle='-',color='gray',label='5sta')
mi30sta=axarr[1].plot(p_mi_sorted_30sta,pmi30pdf,linestyle='--',color='red',label='30sta')
mi32sta=axarr[1].plot(p_mi_sorted_32sta,pmi32pdf,linestyle='-.',color='blue',label='32sta')

axarr[1].hist(p_mi_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='gray',normed=True)
axarr[1].hist(p_mi_sorted_30sta,bins=50,alpha=0.2,histtype='stepfilled',color='red',normed=True)
axarr[1].hist(p_mi_sorted_32sta,bins=50,alpha=0.2,histtype='stepfilled',color='blue',normed=True)

axhandle,axlabel=axarr[0].get_legend_handles_labels()
axarr[1].legend(axhandle,axlabel,loc='upper right',frameon=False)

axarr[1].set_xlabel('Path Term (ln residual')
axarr[1].set_ylabel('Probability Density/Normed count')
axarr[1].title.set_text('Path PDF and Normed histogram for Mixed effects model')

### Save ###
pfig.savefig(home+comp_dir+'/figs/path_pdf.png')


#############################################################3
#############################################################3
#############################################################3
# Make Same figure for stations:

s_tr_sorted_5sta = np.sort(trobj_5sta.site_terms)
s_mi_sorted_5sta = np.sort(miobj_5sta.site_terms)
str5pdf = norm.pdf(s_tr_sorted_5sta)
smi5pdf = norm.pdf(s_mi_sorted_5sta)

s_tr_sorted_30sta = np.sort(trobj_30sta.site_terms)
s_mi_sorted_30sta = np.sort(miobj_30sta.site_terms)
str30pdf = norm.pdf(s_tr_sorted_30sta)
smi30pdf = norm.pdf(s_mi_sorted_30sta)

s_tr_sorted_32sta = np.sort(trobj_32sta.site_terms)
s_mi_sorted_32sta = np.sort(miobj_32sta.site_terms)
str32pdf = norm.pdf(s_tr_sorted_32sta)
smi32pdf = norm.pdf(s_mi_sorted_32sta)

#####################3
#### Plotting ########
######################

# Initiate figure; left for single-mean, right for mixed
fsize=(20,13)
sfig,axarr = plt.subplots(1,2,figsize=fsize)

# First plot the pdf and histogram of the 5sta:
tr5sta=axarr[0].plot(s_tr_sorted_5sta,str5pdf,linestyle='-',color='gray',label='5sta')
tr30sta=axarr[0].plot(s_tr_sorted_30sta,str30pdf,linestyle='--',color='red',label='30sta')
tr32sta=axarr[0].plot(s_tr_sorted_32sta,str32pdf,linestyle='-.',color='blue',label='32sta')

axarr[0].hist(s_tr_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='gray',normed=True)
axarr[0].hist(s_tr_sorted_30sta,bins=50,alpha=0.2,histtype='stepfilled',color='red',normed=True)
axarr[0].hist(s_tr_sorted_32sta,bins=50,alpha=0.2,histtype='stepfilled',color='blue',normed=True)

axhandle,axlabel=axarr[0].get_legend_handles_labels()
axarr[0].legend(axhandle,axlabel,loc='upper right',frameon=False)

axarr[0].set_xlabel('Site Term (ln residual')
axarr[0].set_ylabel('Probability Density/Normed count')
axarr[0].title.set_text('Site PDF and Normed histogram for Single-mean model')

#############
#  Now for mixed:
# First plot the pdf and histogram of the 5sta:
mi5sta=axarr[1].plot(s_mi_sorted_5sta,smi5pdf,linestyle='-',color='gray',label='5sta')
mi30sta=axarr[1].plot(s_mi_sorted_30sta,smi30pdf,linestyle='--',color='red',label='30sta')
mi32sta=axarr[1].plot(s_mi_sorted_32sta,smi32pdf,linestyle='-.',color='blue',label='32sta')

axarr[1].hist(s_mi_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='gray',normed=True)
axarr[1].hist(s_mi_sorted_30sta,bins=50,alpha=0.2,histtype='stepfilled',color='red',normed=True)
axarr[1].hist(s_mi_sorted_32sta,bins=50,alpha=0.2,histtype='stepfilled',color='blue',normed=True)

axhandle,axlabel=axarr[0].get_legend_handles_labels()
axarr[1].legend(axhandle,axlabel,loc='upper right',frameon=False)

axarr[1].set_xlabel('Site Term (ln residual')
axarr[1].set_ylabel('Probability Density/Normed count')
axarr[1].title.set_text('Site PDF and Normed histogram for Mixed effects model')

### Save ###
sfig.savefig(home+comp_dir+'/figs/site_pdf.png')


#############################################################3
#############################################################3
#############################################################3
# Make Same figure for events:

e_tr_sorted_5sta = np.sort(trobj_5sta.E_residual)
e_mi_sorted_5sta = np.sort(miobj_5sta.E_residual)
etr5pdf = norm.pdf(e_tr_sorted_5sta)
emi5pdf = norm.pdf(e_mi_sorted_5sta)

e_tr_sorted_30sta = np.sort(trobj_30sta.E_residual)
e_mi_sorted_30sta = np.sort(miobj_30sta.E_residual)
etr30pdf = norm.pdf(e_tr_sorted_30sta)
emi30pdf = norm.pdf(e_mi_sorted_30sta)

e_tr_sorted_32sta = np.sort(trobj_32sta.E_residual)
e_mi_sorted_32sta = np.sort(miobj_32sta.E_residual)
etr32pdf = norm.pdf(e_tr_sorted_32sta)
emi32pdf = norm.pdf(e_mi_sorted_32sta)

#####################3
#### Plotting ########
######################

# Initiate figure; left for single-mean, right for mixed
fsize=(20,13)
efig,axarr = plt.subplots(1,2,figsize=fsize)

# First plot the pdf and histogram of the 5sta:
tr5sta=axarr[0].plot(e_tr_sorted_5sta,etr5pdf,linestyle='-',color='gray',label='5sta')
tr30sta=axarr[0].plot(e_tr_sorted_30sta,etr30pdf,linestyle='--',color='red',label='30sta')
tr32sta=axarr[0].plot(e_tr_sorted_32sta,etr32pdf,linestyle='-.',color='blue',label='32sta')

axarr[0].hist(e_tr_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='gray',normed=True)
axarr[0].hist(e_tr_sorted_30sta,bins=50,alpha=0.2,histtype='stepfilled',color='red',normed=True)
axarr[0].hist(e_tr_sorted_32sta,bins=50,alpha=0.2,histtype='stepfilled',color='blue',normed=True)

axhandle,axlabel=axarr[0].get_legend_handles_labels()
axarr[0].legend(axhandle,axlabel,loc='upper right',frameon=False)

axarr[0].set_xlabel('Event Term (ln residual')
axarr[0].set_ylabel('Probability Density/Normed count')
axarr[0].title.set_text('Event PDF and Normed histogram for Single-mean model')

#############
#  Now for mixed:
# First plot the pdf and histogram of the 5sta:
mi5sta=axarr[1].plot(e_mi_sorted_5sta,emi5pdf,linestyle='-',color='gray',label='5sta')
mi30sta=axarr[1].plot(e_mi_sorted_30sta,emi30pdf,linestyle='--',color='red',label='30sta')
mi32sta=axarr[1].plot(e_mi_sorted_32sta,emi32pdf,linestyle='-.',color='blue',label='32sta')

axarr[1].hist(e_mi_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='gray',normed=True)
axarr[1].hist(e_mi_sorted_30sta,bins=50,alpha=0.2,histtype='stepfilled',color='red',normed=True)
axarr[1].hist(e_mi_sorted_32sta,bins=50,alpha=0.2,histtype='stepfilled',color='blue',normed=True)

axhandle,axlabel=axarr[0].get_legend_handles_labels()
axarr[1].legend(axhandle,axlabel,loc='upper right',frameon=False)

axarr[1].set_xlabel('Event Term (ln residual')
axarr[1].set_ylabel('Probability Density/Normed count')
axarr[1].title.set_text('Event PDF and Normed histogram for Mixed effects model')

### Save ###
efig.savefig(home+comp_dir+'/figs/event_pdf.png')


#############################################################3
############## 5 sta Single mean and Mixed pdfs ##############
#############################################################3

# Make a single figure (one for event, one for site, one for path) to compare mixed vs. single-mean
#    between just the 5 sta case:

# Events:
ev5stafig=plt.figure()

# pdfs
etr5sta=plt.plot(e_tr_sorted_5sta,etr5pdf,linestyle='--',color='red',label='single-mean')
emi5sta=plt.plot(e_mi_sorted_5sta,emi5pdf,linestyle='-.',color='blue',label='mixed')
# hists
plt.hist(e_tr_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='red',normed=True)
plt.hist(e_mi_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='blue',normed=True)

# Get z test value???
e5sta_z = ra.z_test(e_tr_sorted_5sta,e_mi_sorted_5sta)


plt.legend(loc='upper right',frameon=False)

plt.xlabel('Event Term (ln residual)')
plt.ylabel('Probability/Normed count')
plt.title('Event Terms: Single-mean vs. Mixed effects \n Z-test: %f'  % e5sta_z)

ev5stafig.savefig(home+comp_dir+'/figs/event_5sta.png')

######
# Stations:
st5stafig=plt.figure()

# pdfs
str5sta=plt.plot(s_tr_sorted_5sta,str5pdf,linestyle='--',color='red',label='single-mean')
smi5sta=plt.plot(s_mi_sorted_5sta,smi5pdf,linestyle='-.',color='blue',label='mixed')
# hists
plt.hist(s_tr_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='red',normed=True)
plt.hist(s_mi_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='blue',normed=True)

# Get z test value???
s5sta_z  = ra.z_test(s_tr_sorted_5sta,s_mi_sorted_5sta)


plt.legend(loc='upper right',frameon=False)

plt.xlabel('Site Term (ln residual)')
plt.ylabel('Probability/Normed count')
plt.title('Site Terms: Single-mean vs. Mixed effects \n Z-test: %f' % s5sta_z)

st5stafig.savefig(home+comp_dir+'/figs/site_5sta.png')


######
# Paths:
pa5stafig=plt.figure()

# pdfs
ptr5sta=plt.plot(p_tr_sorted_5sta,ptr5pdf,linestyle='--',color='red',label='single-mean')
pmi5sta=plt.plot(p_mi_sorted_5sta,pmi5pdf,linestyle='-.',color='blue',label='mixed')
# hists
plt.hist(p_tr_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='red',normed=True)
plt.hist(p_mi_sorted_5sta,bins=50,alpha=0.2,histtype='stepfilled',color='blue',normed=True)

# Get z test value???
p5sta_z = ra.z_test(p_tr_sorted_5sta,p_mi_sorted_5sta)

plt.legend(loc='upper right',frameon=False)

plt.xlabel('Path Term (ln residual)')
plt.ylabel('Probability/Normed count')
plt.title('Path Terms: Single-mean vs. Mixed effects \n Z-test: %f' % p5sta_z)

pa5stafig.savefig(home+comp_dir+'/figs/path_5sta.png')



#############################################################3
#############################################################3
#############################################################3
#z_test_5_32_tr_mi_path=(np.mean(trobj_5sta.path_terms) - mean(miobj_5sta.path_terms))/np.sqrt(std(trobj_5sta.path_terms)**2 + std(miobj_5sta.path_terms)**2)