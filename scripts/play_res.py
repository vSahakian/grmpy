import dread as dr
import cdefs as cdf
import numpy as np
import inversion as inv
import pickle
import matplotlib.pyplot as plt
import gmpe as gm
import rescomp as rcomp


#Change this parameter depending on where you run:
#0=desktop
#1=mac

what_home=0

if what_home==0:
    #Desktop:
    HOME='/home/vsahakian'
elif what_home==1:
    #Mac:
    HOME='/Users/vsahakian'



#AB's flatfile:
ffile=HOME+'/anza/data/Anzadata_Acc_Vel_May2016_40_50.mat'

#Location to store figures:
fig_dir=HOME+'/anza/models/figs/'
obj_dir=HOME+'/anza/models/pckl/'

#Read data from .mat file, store important values:
ev,sta,N,ml,mw,DA,DV,r,vs30=dr.mread(ffile)

#Load into database object
abdb=cdf.db(ev,sta,N,ml,mw,DA,DV,r,vs30)


##Residual Computation####
#Load in the model to use:
fname=obj_dir+'regr_0.0_6.5_resid_2676.06963031.pckl'
model=pickle.load(open(fname,'r'))

#Overall residual, 
#In some places, vs30 is 0.  Set these to vref.
vref=760
mdep_ffdf=0
#Where are they 0?
vs30_0ind=np.where(abdb.vs30==0)[0]
#Get vs30 from database...
vs30=abdb.vs30
#Set 0 entries to vref:
vs30[vs30_0ind]=vref

#Now compute the predicted value of PGA...
d_predicted=gm.compute_model(model.m,model.rng,abdb.mw,abdb.r,abdb.ffdf,vs30,vref,mdep_ffdf)

#Get residuals:
total_residuals,mean_residual,std_dev=rcomp.total_residual(abdb,d_predicted)


#Get unique events:
unique_events=np.unique(abdb.evnum)

#Make a class for each event; append them to a list, made empty here:
event_list=[]
d_predicted_list=[]

#Loop through the unique events, make each into an object, append to event list
for i in range(len(unique_events)):
    unique_ind=np.where(unique_events[i]==abdb.evnum)[0]
    #Get the predicted data for this event only, for each recording
    #d_predicted_i is an array, with len = # of recordings for event):
    d_predicted_i=d_predicted[unique_ind]
    
    #Get the database info for this event:
    evnum_i=abdb.evnum[unique_ind]
    sta_i=abdb.sta[unique_ind]
    stnum_i=abdb.stnum[unique_ind]
    ml_i=abdb.ml[unique_ind]
    mw_i=abdb.mw[unique_ind]
    pga_i=abdb.pga[unique_ind]
    pgv_i=abdb.pgv[unique_ind]
    pga_pg_i=abdb.pga_pg[unique_ind]
    r_i=abdb.r[unique_ind]
    ffdf_i=abdb.ffdf[unique_ind]
    md_ffdf_i=abdb.md_ffdf[unique_ind]
    
    #The only one that is not directly from teh database is vs30; this is because
    #there were some 0's in it, so above that is compensated for by changing 
    #them to vref.
    vs30_i=vs30[unique_ind]
    
    eventi=cdf.event(evnum_i,sta_i,stnum_i,ml_i,mw_i,pga_i,pgv_i,pga_pg_i,r_i,vs30_i,ffdf_i,md_ffdf_i)
    
    #Apped the event object, and the d_predicted to the list:
    event_list.append(eventi)
    d_predicted_list.append(d_predicted_i)
 
    
 ###Get Event Residuals###
    
#Then for each event, compute the total residual.
#Zero out arrays...
E_evnum=[]
E_mw=[]
E_residual=[]
E_std_dev=[]

#Loop over each event object:
for eventi in range(len(event_list)):

    #Get the event residual for each event:
    evnum_i,evmw_i,E_residual_i,std_dev_i=rcomp.event_residual(event_list[eventi],d_predicted_list[eventi])
    
    #Append to the residual arrays:
    E_evnum.append(evnum_i)
    E_mw.append(evmw_i)
    E_residual.append(E_residual_i)
    E_std_dev.append(std_dev_i)
    
E_evnum=np.array(E_evnum)
E_mw=np.array(E_mw)
E_residual=np.array(E_residual)
E_std_dev=np.array(E_std_dev)
    

###Get Within-Event Residuals###

#Zero out the arrays/lists:
W_evnum=[]
W_mw=[]
W_sta=[]
W_stnum=[]
W_residuals=[]
W_mean=[]
W_std_dev=[]

#Loop over each event object, again - this time provide event residual as well:
for eventi in range(len(event_list)):
    
    #Get the Within-Event residual and other info for each event:
    evnum_i,evmw_i,sta_i,stnum_i,W_residuals_i,W_mean_i,W_std_dev_i=rcomp.within_event_residual(event_list[eventi],d_predicted[eventi],E_residual[eventi])
    
    #Append to the lists/arrays:
    W_evnum.append(evnum_i)
    W_mw.append(evmw_i)
    W_sta.append(sta_i)
    W_stnum.append(stnum_i)
    W_residuals.append(W_residuals_i)
    W_mean.append(W_mean_i)
    W_std_dev.append(W_std_dev_i)
    
#Save to arrays:
W_evnum=np.array(W_evnum)
W_mw=np.array(W_mw)
#Don't convert W_sta...it's a bunch of strings...
W_stnum=np.array(W_stnum)
W_residuals=np.array(W_residuals) 
W_mean=np.array(W_mean)
W_std_dev=np.array(W_std_dev)   
    
    