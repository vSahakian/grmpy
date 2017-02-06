#Testing file....

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
invrun='test_plotting'
fig_dir=home+'/models/figs/'+invrun+'/'
obj_dir=home+'/models/pckl/'+invrun+'/'
model_dir=home+'/models/pckl/'+invrun+'/'


#Filename:
dbpath=home+'/data/databases/v2anza2013/v2anza2013_pgrid_5sta_res4.pckl'

print 'Using database %s' % dbpath

####################################################
############Inversion Parameters####################
####################################################

#Coefficient file for ASK
coeff_file=coeff_file=home+'/data/coeffs/ASK14_coeffs.m'

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
bmax=120

#Number of coefficients to solve for in functional form:
ncoeff=5

#Axis limits for plotting:
axlims=[[0,6],[-8,-0.5]]

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
plotdist=np.array([0,10,20,40,90,150,200])

# plotting distance for ASK:
ask_dist=10


#####################################################
###########Run Inversion and plot:###################
#####################################################
#
#
#################
######Setup######
#################

#Get the string for the filename, based on the ranges:
for k in range(len(rng)):
    if k==0:
        strname=np.str(rng[k])
    else:
        strname=strname+'_'+np.str(rng[k])

########
#Invert#
########

print 'Setting up inversion'
inv_dat=run_inv.setup_run_inversion(home,dbpath,dbname,ncoeff,rng,sdist,Mc,smth,vref,mdep_ffdf)

print inv_dat

#Get basename for model:
basename='regr_Mc'+str(Mc)+'_'+strname+'_VR_'+str(np.around(inv_dat.VR,decimals=1))
#basename='regr_'+strname+'_VR_'+np.str(np.around(inv_dat.VR,decimals=1))
modelpath=model_dir+basename+'.pckl'

print 'Will read in '+modelpath

print 'plotting...'
#Plot:
fig1=run_inv.plot_data_model(home,dbpath,dbname,modelpath,coeff_file,mdep_ffdf,plotdist,ask_dist,Mc,axlims,bmin,bmax,vref)



###################################################################################
#################                 MIXED EFFECTS                 ###################
###################################################################################
#
#print 'Running mixed effects'
##Now try with mixed effects:
#run_name = 'mixedregr_v2anza2013_Mc_8.5_res4'
#run_home=home+'/models/residuals/'
#resaxlim_r = [[0,180],[-5,5]]
#resaxlim_mw = [[0,4],[-5,5]]
#
##Fictitious depth parameter:
#c=4.5
#
##Initialize the residuals directories:
#inithome=HOME+'/anza/models/residuals/'
#
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
#    
#    
## Now run mixed effects approach #
#invdat,invpath,tresid,mixed_residuals,d_r_prediction,mixed_resid_path=run_inv.run_mixedeffects(home,codehome,run_name,dbpath,dbname,Mc,vref,c)
#
#print 'Plotting mixed effects model with data'
## Plot data with model:
#mixedinv = run_inv.plot_data_model(home,dbpath,dbname,invpath,coeff_file,mdep_ffdf,plotdist,Mc,axlims,bmin,bmax,vref)
#
## Get some mean values for the stats file:
#mean_tot=np.mean(mixed_residuals.total_residual)
#std_dev_tot=np.std(mixed_residuals.total_residual)
#
#E_mean=np.mean(mixed_residuals.E_residual)
#E_std_dev=np.std(mixed_residuals.E_residual)
#
#print 'plotting all residuals'
## Plot all residuals:
#run_res.plot_total(tresid,home,run_name,resaxlim_mw)
#
#print 'plotting event terms'
## Plot event terms:
#eventfig1 = run_res.makeEvents_mixed(run_home,run_name,mixed_resid_path,Mc,vref,mdep_ffdf,resaxlim_mw)
#
## Get station objects:
#station_list=run_res.makeStations_mixed(run_home,run_name,mixed_resid_path)
#
## Plot W residuals on one plot:
#W_mean,W_std_dev = run_res.plot_Wresid(run_home,run_name,resaxlim_mw)
#
#print 'plotting station terms'
## Plot by station:
#run_res.plot_site_WE(run_home,run_name,resaxlim_mw)
#
#print 'plotting path residuals'
## Plot path residuals:
#f_mw,f_dist,pterm_mean,pterm_std=run_res.plotPath_mixed(run_home,run_name,mixed_resid_path,resaxlim_mw,resaxlim_r)
#
#print 'Writing mixed stats to a file'
## Writing stats to a file
#run_res.write_stats_mixed(run_home,run_name,mixed_resid_path,mean_tot,std_dev_tot,E_mean,E_std_dev,W_mean,W_std_dev,pterm_mean,pterm_std)