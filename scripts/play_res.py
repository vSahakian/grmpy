import dread as dr
import cdefs as cdf
import numpy as np
import inversion as inv
import pickle
import matplotlib.pyplot as plt
import gmpe as gm


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


#Get unique events:
unique_events=np.unique(abdb.evnum)

#Make a class for each event; append them to a list, made empty here:
event_list=[]

#Loop through the unique events, make each into an object, append to event list
for i in range(len(unique_events)):
    unique_ind=np.where(unique_events[i]==abdb.evnum)[0]
    eventi=cdf.event(abdb.evnum,abdb.sta,abdb.stnum,abdb.ml,abdb.mw,abdb.pga,abdb.pgv,abdb.pga_pg,abdb.r,abdb.vs30,abdb.ffdf,abdb.md_ffdf)
    event_list.append(eventi)
    
#Then for each event, compute the total residual:
for eventi in range(len(event_list)):
    
