#Testing file....

import dread as dr
import cdefs as cdf
import numpy as np
import inversion as inv

ffile='/Users/vsahakian/anza/data/Anzadata_Acc_Vel_May2016_40_50.mat'

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


#Inverting:
rng=np.array([0,3.3,6.5])
sdist=np.array([1,5,10,15,20])
smth=100

G,d=inv.iinit_pga(abdb,rng,sdist,smth)