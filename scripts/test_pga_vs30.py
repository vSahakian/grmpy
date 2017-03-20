# Test the PGV and ln inversion of data...

import numpy as np
import cPickle as pickle
import run_inversion as run_inv
import run_res


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

#Location to store figures:
home=HOME+'/anza'
invrun='v2anza2013_pga_vs30'
fig_dir=home+'/models/figs/'+invrun+'/'
obj_dir=home+'/models/pckl/'+invrun+'/'
model_dir=home+'/models/pckl/'+invrun+'/'

#Coefficient file for ASK
coeff_file=coeff_file=home+'/data/coeffs/ASK14_coeffs.m'

#Filename:
dbpath=home+'/data/databases/v2anza2013/v2anza2013_pgrid_5sta_res4.pckl'

print 'Using database %s' % dbpath

####################################################
############Inversion Parameters####################
####################################################

#Path name:
dbname=invrun

# Invert for?
predictive_parameter='pga'

# Correct the data by vs30?
vs30_correct=0    # 0 = no correction, 1 = yes, vs30 correction


#Define the ranges for the inversion - at each range boundary (i.e., between
#[0:3.3], and [3.3:4.5], the solution will be smoothed so there is no jump at 
#the range boundary
rng=np.array([0,6.5])

#Reference Vs30:
vref=760

#Bin min and max for data scatter colorbar (distances):
bmin=10
bmax=120

#Number of coefficients to solve for in functional form:
ncoeff=6

#Axis limits for plotting:
#axlims=[[0,6],[-8,-0.5]]
axlims=[[0,5.5],[-7,-0.5]]


#Number of distances to include in smoothing - there will be this many extra
#equations added on at each range boundary
smin=0
smax=280
sstep=1
sdist=np.array(range(smin,smax,sstep))

##Centering Magnitude for x**2 term
Mc=8.5

#Smoothing factor
smth=500

#Use magnitude-dependent fictitious depth/finite fault dimension factor?####
#no == 0, yes == 1
mdep_ffdf=0

#Plotting distances:
#plotdist=np.array([0,20,40,60,120])
plotdist=np.array([0,10,20,40,80,160,220])

# plotting distance for ASK:
ask_dist=10

## For plotting only:
setmodel='/Users/vsahakian/anza/models/pckl/v2anza2013_pga_vs30/regr_pga_Mc8.5_0.0_6.5_VR_99.4.pckl'
#setmixedmodel='/Users/vsahakian/anza/models/pckl/v2anza2013/mixedregr_v2anza2013_Mc_8.5_VR_99.9.pckl'


#################
######Setup######
#################

#Get the string for the filename, based on the ranges:
for k in range(len(rng)):
    if k==0:
        strname=np.str(rng[k])
    else:
        strname=strname+'_'+np.str(rng[k])

#########
##Invert#
#########
#
#print 'Setting up inversion'
#inv_dat=run_inv.setup_run_inversion(home,dbpath,dbname,ncoeff,rng,sdist,Mc,smth,vref,mdep_ffdf,predictive_parameter=predictive_parameter,data_correct=vs30_correct)
#
#
## Plot the data and model now:
#
#
## If you just want to plot with a certain model:
#modelpath=setmodel
#fig1=run_inv.plot_data_model(home,dbpath,dbname,modelpath,coeff_file,mdep_ffdf,plotdist,ask_dist,Mc,axlims,bmin,bmax,vref,predictive_parameter=predictive_parameter,ncoeff=ncoeff,data_correct=0)



###################################################################################
#################                 MIXED EFFECTS                 ###################
###################################################################################
#
print 'Running mixed effects'
#Now try with mixed effects:
run_name = 'mixedregr_pga_v2anza2013_Mc_8.5_res4'
run_home=home+'/models/residuals/'
resaxlim_r = [[0,180],[-5,5]]
resaxlim_mw = [[0,4],[-5,5]]


#Fictitious depth parameter:
c=4.5

#Initialize the residuals directories:
inithome=HOME+'/anza/models/residuals/'

#runall=1
#
##Initialize directories:
#runall=run_res.init(inithome,run_name)
#
#if runall==0:
#    print 'Not clobbering, exiting...'
#    
#elif runall==1:
#    print 'Continuing...'
    
    
# Now run mixed effects approach #
invdat,invpath,tresid,mixed_residuals,d_r_prediction,mixed_resid_path=run_inv.run_mixedeffects(home,codehome,run_name,dbpath,dbname,Mc,vref,c,predictive_parameter='pga',ncoeff=5,data_correct=0,a1='none',a2='none',a3='none',a4='none',a5='none',a6='none')


