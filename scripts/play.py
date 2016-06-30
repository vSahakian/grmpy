#Testing file....

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

#Location to store figures:
fig_dir=HOME+'/anza/models/figs/'
obj_dir=HOME+'/anza/models/pckl/'

######Making the database object - once made and saved, do not run this again######
#
##AB's flatfile:
#ffile=HOME+'/anza/data/Anzadata_Acc_Vel_May2016_40_50.mat'
#hfile=HOME+'/anza/data/YSH_2010.hash'
#
##Read data from .mat file, store important values:
#ev,sta,N,ml,mw,DA,DV,r,vs30,lat,lon,depth=dr.mread(ffile,hfile)
#
##Load into database object
#abdb=cdf.db(ev,sta,N,ml,mw,DA,DV,r,vs30,lat,lon,depth)
#
##Save the database object:
#fname=HOME+'/anza/data/abdb.pckl'
#datobj=open(fname,'w')
#pickle.dump(abdb,datobj)
#datobj.close()
########################End making database object############################


##########Open the database object:###################
#Filename:
fname=HOME+'/anza/data/abdb.pckl'
datobj=open(fname,'r')
abdb=pickle.load(datobj)
datobj.close()

#Plot pga's:
#abdb.plot_apga()

#Plot by distance:
bmin=0
bmax=22
step=1
#abdb.plot_rpga(bmin,bmax,step)


#Number of coefficients to solve for in functional form:
ncoeff=5
#Define the ranges for the inversion - at each range boundary (i.e., between
#[0:3.3], and [3.3:4.5], the solution will be smoothed so there is no jump at 
#the range boundary
rng=np.array([0,3.3,4.5,6.5])
#rng=np.array([0,1,1.5,2,2.5,3.3,6.5])
#rng=np.array([0,1,1.5,2,2.5,3.5,6.5])
#rng=np.array([0,1.5,3,4,6.5])
rng=np.array([0,6.5])#####Making the database object - once made and saved, do not run this again######

#AB's flatfile:
ffile=HOME+'/anza/data/Anzadata_Acc_Vel_May2016_40_50.mat'
hfile=HOME+'/anza/data/YSH_2010.hash'

#Location to store figures:
fig_dir=HOME+'/anza/models/figs/'
obj_dir=HOME+'/anza/models/pckl/'

#Read data from .mat file, store important values:
ev,sta,N,ml,mw,DA,DV,r,vs30,lat,lon,depth=dr.mread(ffile,hfile)

#Load into database object
abdb=cdf.db(ev,sta,N,ml,mw,DA,DV,r,vs30,lat,lon,depth)

#Save the database object:
fname=HOME+'/anza/data/abdb.pckl'
datobj=open(fname,'w')
pickle.dump(abdb,datobj)
datobj.close
#rng=np.array([0,1,2,3.3,4.5,6.5])
#rng=np.array([0,6.5])
#rng=np.array([0,2,3,4,6.5])


#Number of distances to include in smoothing - there will be this many extra
#equations added on at each range boundary
sdist=np.array([1,5,10,15,20])

#Smoothing factor
smth=500

#Use magnitude-dependent fictitious depth/finite fault dimension factor?####
#no == 0, yes == 1
mdep_ffdf=0

#Invert:
#Make matrices
G,d=inv.iinit_pga(abdb,ncoeff,rng,sdist,smth,mdep_ffdf)
#Invert
m, resid, L2norm, VR, rank, svals=inv.invert(G,d)

#Compute the predicted value (from the GMPE) at each data point
vref=760
#Magnitude dependent fictitous depth?
if mdep_ffdf==0:
    d_predicted=gm.compute_model(m,rng,abdb.mw,abdb.r,abdb.ffdf,abdb.vs30,vref,mdep_ffdf)
elif mdep_ffdf==1:
    d_predicted=gm.compute_model(m,rng,abdb.mw,abdb.r,abdb.md_ffdf,abdb.vs30,vref,mdep_ffdf)
    

#Compute the magnitude/log10pga for each distance, to plot on top of data:
mw_model,d_model=gm.compute_model_fixeddist(m,rng,sdist,mdep_ffdf)

#Get the NGA predictions to plot on the same figure:
#Coefficient file:
coeff_file=coeff_file=HOME+'/anza/data/coeffs/ASK14_coeffs.m'
#Do it just for one distance for now, say R=5km.  
Rrup=5*np.ones(abdb.r.shape)
#Get the NGA predictions...
f1,M_sort,f1_sort=gm.ask2014_pga(abdb.mw,Rrup,coeff_file,1,[0,0])



#Plotting params...
axlims=[[1,6],[-7,0]]
#Plot against data to check:
abdb.plot_rpga_withmodel(bmin,bmax,step,mw_model,d_model,rng,sdist,axlims,VR,vref)
#Plot NGA:
plt.plot(M_sort,f1_sort,'--',label='ASK2014')


#Save plots:
#Get the string for the filename, based on the ranges:

for k in range(len(rng)):
    if k==0:
        strname=np.str(rng[k])
    else:
        strname=strname+'_'+np.str(rng[k])
    
basename='regr_'+strname+'_VR_'+np.str(np.around(VR,decimals=1))
figname=fig_dir+basename+'.png'
plt.savefig(figname,transparent=True)

#Save G, d, and m.....and other things...
#Put into an inversion object:
invdat=cdf.invinfo(G,d,m,resid,L2norm,VR,rank,svals,rng,sdist,smth)
fname=obj_dir+basename+'.pckl'
datobj=open(fname,'w')
pickle.dump(invdat,datobj)
datobj.close()

