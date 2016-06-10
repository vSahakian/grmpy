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
#rng=np.array([0,3.3,6.5])
#sdist=np.array([1,5,10,15,20])
#smth=100

#Number of coefficients to solve for in functional form:
ncoeff=5
#Define the ranges for the inversion - at each range boundary (i.e., between
#[0:3.3], and [3.3:4.5], the solution will be smoothed so there is no jump at 
#the range boundary
rng=np.array([0,3.3,4.5,6.5])
#Number of distances to include in smoothing - there will be this many extra
#equations added on at each range boundary
sdist=np.array([1,5,10,15,20])
#Smoothing factor
smth=100

G,d=inv.iinit_pga(abdb,ncoeff,rng,sdist,smth)
m, resid, rank, svals=inv.invert(G,d)

#Plot against data to check:
abdb.plot_rpga_withmodel(bmin,bmax,step,m,rng,sdist)

np.save(
