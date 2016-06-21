#Testing file....

import dread as dr
import cdefs as cdf
import numpy as np
import inversion as inv
import pickle
import matplotlib.pyplot as plt


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
#rng=np.array([0,6.5])
#rng=np.array([0,1,2,3.3,4.5,6.5])
#rng=np.array([0,6.5])



#Number of distances to include in smoothing - there will be this many extra
#equations added on at each range boundary
sdist=np.array([1,5,10,15,20])

#Smoothing factor
smth=500

#Use magnitude-dependent fictitious depth/finite fault dimension factor?
#no == 0, yes == 1
mdep_ffdf=1

#Invert:
#Make matrices
G,d=inv.iinit_pga(abdb,ncoeff,rng,sdist,smth,mdep_ffdf)
#Invert
m, resid, L2norm, VR, rank, svals=inv.invert(G,d)



#m=np.array([  2.39195267e+00,   1.33068285e+00,   2.20961871e-18,
#        9.00000000e-01,  -1.06447647e-01,   1.37851521e+01,
#       -9.26696628e-01,  -1.45862809e-01,   9.00000000e-01,
#       -1.06430279e-01,   6.38186568e+00,   5.41230552e-01,
#       -9.60255039e-02,   9.00119541e-01,  -1.06435188e-01])

#Plotting params...
vref=760
axlims=[[1,6],[-7,0]]
#Plot against data to check:
abdb.plot_rpga_withmodel(bmin,bmax,step,m,rng,sdist,axlims,VR,vref)

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

