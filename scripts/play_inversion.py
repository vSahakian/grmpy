#Testing file....

import numpy as np
import cPickle as pickle
import run_inversion as run_inv


#Change this parameter depending on where you run:
#0=desktop
#1=mac

what_home=1

if what_home==0:
    #Desktop:
    HOME='/media/vsahakian/katmai'
elif what_home==1:
    #Mac:
    HOME='/Users/vsahakian'

#Location to store figures:
home=HOME+'/anza'
fig_dir=home+'/models/figs/'
obj_dir=home+'/models/pckl/'

#AB's flatfile:
ffile=home+'/data/Anzadata_Acc_Vel_May2016_40_50.mat'
hfile=home+'/data/YSH_2010.hash'

#Coefficient file for ASK
coeff_file=coeff_file=home+'/data/coeffs/ASK14_coeffs.m'


###################################################################################
######Making the database object - once made and saved, do not run this again######
###################################################################################

#
##AB's flatfile:
#ffile=HOME+'/anza/data/Anzadata_Acc_Vel_May2016_40_50.mat'
#hfile=HOME+'/anza/data/YSH_2010.hash'
#sfile=HOME+'/anza/data/stations/AZ.ll'
#
##Read data from .mat file, store important values:
#ev,sta,N,ml,mw,DA,DV,r,vs30,lat,lon,depth,stlat,stlon,source_i,receiver_i=dr.mread(ffile,hfile,sfile)
#
##Load into database object
#abdb=cdf.db(ev,sta,N,ml,mw,DA,DV,r,vs30,lat,lon,depth,stlat,stlon,source_i,receiver_i)
#
##Save the database object:
#fname=HOME+'/anza/data/abdb.pckl'
#datobj=open(fname,'w')
#pickle.dump(abdb,datobj)
#datobj.close()

##############################################################################
########################End making database object############################
##############################################################################


####################################################
##########Run Inversion and plot:###################
####################################################

#Open the database object;

#Filename:
##
#fname=HOME+'/anza/data/abdb.pckl'
#datobj=open(fname,'r')
#abdb=pickle.load(datobj)
#datobj.close()
##

#Filename:
dbpath=home+'/data/abdb_5sta.pckl'
datobj=open(dbpath,'r')
abdb=pickle.load(datobj)
datobj.close()

#Define the ranges for the inversion - at each range boundary (i.e., between
#[0:3.3], and [3.3:4.5], the solution will be smoothed so there is no jump at 
#the range boundary
#rng=np.array([0,3.3,4.5,6.5])
#rng=np.array([0,1,1.5,2,2.5,3.3,6.5])
#rng=np.array([0,1,1.5,2,2.5,3.5,6.5])
#rng=np.array([0,1.5,3,4,6.5])
#rng=np.array([0,3.3,6.5])
#rng=np.array([0,3,6.5])
#rng=np.array([0,1,2,6.5])
rng=np.array([0,6.5])


#####################
######Inversion######
#####################

#Reference Vs30:
vref=760

#Bin min and max for data scatter colorbar:
bmin=0
bmax=22

#Number of coefficients to solve for in functional form:
ncoeff=5

#Axis limits for plotting:
axlims=[[1,6],[-7,0]]

#Number of distances to include in smoothing - there will be this many extra
#equations added on at each range boundary
sdist=np.array([1,5,10,15,20])

#Smoothing factor
smth=500

#Use magnitude-dependent fictitious depth/finite fault dimension factor?####
#no == 0, yes == 1
mdep_ffdf=0

#Get the string for the filename, based on the ranges:
for k in range(len(rng)):
    if k==0:
        strname=np.str(rng[k])
    else:
        strname=strname+'_'+np.str(rng[k])

########
#Invert#
########

inv_dat=run_inv.setup_run_inversion(home,dbpath,ncoeff,rng,sdist,smth,mdep_ffdf)

print inv_dat

#Get basename for model:
basename='regr_'+strname+'_VR_'+np.str(np.around(inv_dat.VR,decimals=1))
modelpath=home+'/models/pckl/'+basename+'.pckl'

#Plot:
fig1=run_inv.plot_data_model(home,dbpath,modelpath,coeff_file,mdep_ffdf,sdist,axlims,bmin,bmax,vref)




