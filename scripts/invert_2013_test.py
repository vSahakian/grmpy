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
invrun='test2013/'
fig_dir=home+'/models/figs/'+invrun
obj_dir=home+'/models/pckl/'+invrun
model_dir=home+'/models/pckl/'+invrun

#Coefficient file for ASK
coeff_file=coeff_file=home+'/data/coeffs/ASK14_coeffs.m'

#Filename:
dbpath=home+'/data/databases/db2013_test/db2013test_5sta.pckl'
#Path name:
dbname=invrun

#Define the ranges for the inversion - at each range boundary (i.e., between
#[0:3.3], and [3.3:4.5], the solution will be smoothed so there is no jump at 
#the range boundary
rng=np.array([0,6.5])

#Reference Vs30:
vref=760

#Bin min and max for data scatter colorbar (distances):
bmin=10
bmax=240

#Number of coefficients to solve for in functional form:
ncoeff=5

#Axis limits for plotting:
axlims=[[0,6],[-10,-2]]

#Number of distances to include in smoothing - there will be this many extra
#equations added on at each range boundary
smin=0
smax=280
sstep=10
sdist=np.array(range(smin,smax,sstep))

#Smoothing factor
smth=500

#Use magnitude-dependent fictitious depth/finite fault dimension factor?####
#no == 0, yes == 1
mdep_ffdf=0

#Plotting distances:
#plot_min=0
#plot_max=300
#plot_step=60
#plotdist=np.array(range(plot_min,plot_max,plot_step))
plotdist=np.array([0,20,40,60,160,300])


####################################################
##########Run Inversion and plot:###################
####################################################

##Open the database object;
#datobj=open(dbpath,'r')
#db=pickle.load(datobj)
#datobj.close()


#####################
######Inversion######
#####################



#Get the string for the filename, based on the ranges:
for k in range(len(rng)):
    if k==0:
        strname=np.str(rng[k])
    else:
        strname=strname+'_'+np.str(rng[k])

########
#Invert#
########

inv_dat=run_inv.setup_run_inversion(home,dbpath,dbname,ncoeff,rng,sdist,smth,mdep_ffdf)

print inv_dat

#Get basename for model:
basename='regr_'+strname+'_VR_'+np.str(np.around(inv_dat.VR,decimals=1))
modelpath=model_dir+basename+'.pckl'

#Plot:
fig1=run_inv.plot_data_model(home,dbpath,dbname,modelpath,coeff_file,mdep_ffdf,plotdist,axlims,bmin,bmax,vref)




