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
invrun='v3anza2013_pga_noVs30_5coeff'
fig_dir=home+'/models/figs/'+invrun+'/'
obj_dir=home+'/models/pckl/'+invrun+'/'
model_dir=home+'/models/pckl/'+invrun+'/'

#Coefficient file for ASK
coeff_file=coeff_file=home+'/data/coeffs/ASK14_coeffs.m'

#Filename:
dbpath=home+'/data/databases/v3anza2013/v3anza2013_pgrid_5sta_res4.pckl'

print 'Using database %s' % dbpath

####################################################
############Inversion Parameters####################
####################################################

#Path name:
dbname=invrun

# Invert for?
predictive_parameter='pga'

# Correct the data by vs30?
vs30_correct=0    # 0 = no correction, !0 = yes, vs30 correction


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
#setmodel='/Users/vsahakian/anza/models/pckl/v2anza2013_pga_vs30_%icoeff/regr_pga_Mc8.5_0.0_6.5_VR_99.4.pckl' % ncoeff
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
inv_dat=run_inv.setup_run_inversion(home,dbpath,dbname,ncoeff,rng,sdist,Mc,smth,vref,mdep_ffdf,predictive_parameter=predictive_parameter,data_correct=vs30_correct)

print 'Variance reduction is %f' % inv_dat.VR

#Get basename for model:
basename='regr_' + predictive_parameter + '_Mc'+str(Mc)+'_'+strname+'_VR_'+str(np.around(inv_dat.VR,decimals=1))
#basename='regr_'+strname+'_VR_'+np.str(np.around(inv_dat.VR,decimals=1))
modelpath=model_dir+basename+'.pckl'

# Plot the data and model now:
fig1=run_inv.plot_data_model(home,dbpath,dbname,modelpath,coeff_file,mdep_ffdf,plotdist,ask_dist,Mc,axlims,bmin,bmax,vref)

## If you just want to plot with a certain model:
#modelpath=setmodel
#fig1=run_inv.plot_data_model(home,dbpath,dbname,modelpath,coeff_file,mdep_ffdf,plotdist,ask_dist,Mc,axlims,bmin,bmax,vref,predictive_parameter=predictive_parameter,ncoeff=ncoeff,data_correct=0)



###################################################################################
#################         MIXED EFFECTS invert for all          ###################
###################################################################################
#
print 'Running mixed effects'
#Now try with mixed effects:
run_name = 'mixedregr_v3anza2013_pga_%icoeff_Mc_8.5_res4_noVs30' % ncoeff
run_home=home+'/models/residuals/'
resaxlim_r = [[0,180],[-5,5]]
resaxlim_mw = [[0,4],[-5,5]]


#Fictitious depth parameter:
c=4.5

#Initialize the residuals directories:
inithome=HOME+'/anza/models/residuals/'

runall=1

#Initialize directories:
runall=run_res.init(inithome,run_name)

if runall==0:
    print 'Not clobbering, exiting...'
    
elif runall==1:
    print 'Continuing...'
    
    
    # Now run mixed effects approach #
    invdat,invpath,tresid,mixed_residuals,d_r_prediction,mixed_resid_path=run_inv.run_mixedeffects(home,codehome,run_name,dbpath,dbname,Mc,vref,c,predictive_parameter='pga',ncoeff=ncoeff,data_correct=vs30_correct,a1='none',a2='none',a3='none',a4='none',a5='none',a6='none')
    
    print 'Variance reduction is %f' % invdat.VR
    
    print 'Plotting mixed effects model with data'
    # Plot data with model:
    mixedinv = run_inv.plot_data_model(home,dbpath,dbname,invpath,coeff_file,mdep_ffdf,plotdist,ask_dist,Mc,axlims,bmin,bmax,vref,predictive_parameter=predictive_parameter,ncoeff=ncoeff,data_correct=vs30_correct)

    ## If just plotting:
    #invpath=setmixedmodel
    #mixedinv = run_inv.plot_data_model(home,dbpath,dbname,invpath,coeff_file,mdep_ffdf,plotdist,ask_dist,Mc,axlims,bmin,bmax,vref)


    # Get some mean values for the stats file:
    mean_tot=np.mean(mixed_residuals.total_residual)
    std_dev_tot=np.std(mixed_residuals.total_residual)
    
    # Get the event mean and stadndard deviation
    E_mean=mixed_residuals.E_mean
    E_std_dev=mixed_residuals.E_std

    print 'plotting all residuals'
    # Plot all residuals:
    run_res.plot_total(tresid,home,run_name,resaxlim_mw)

    print 'plotting event terms'
    # Plot event terms:
    eventfig1 = run_res.makeEvents_mixed(run_home,run_name,mixed_resid_path,Mc,vref,mdep_ffdf,resaxlim_mw)

    # Get station objects:
    station_list=run_res.makeStations_mixed(run_home,run_name,mixed_resid_path)

    # Plot W residuals on one plot:
    W_mean,W_std_dev = run_res.plot_Wresid(run_home,run_name,resaxlim_mw)

    print 'plotting station terms'
    # Plot by station:
    run_res.plot_site_WE(run_home,run_name,resaxlim_mw)

    print 'plotting path residuals'
    # Plot path residuals:
    f_mw,f_dist,pterm_mean,pterm_std=run_res.plotPath_mixed(run_home,run_name,mixed_resid_path,resaxlim_mw,resaxlim_r)
    
    print 'Writing mixed stats to a file'
    # Writing stats to a file
    run_res.write_stats_mixed(run_home,run_name,mixed_resid_path,mean_tot,std_dev_tot,E_mean,E_std_dev,W_mean,W_std_dev,pterm_mean,pterm_std)
    
    
    
    
###################################################################################
#################       MIXED EFFECTS with fixed a2 - a5        ###################
###################################################################################
#

## Read in model so you don't have to run the entire thing:
#fixedmodel = '/Users/vsahakian/anza/models/pckl/v3anza2013_pga_vs30_5coeff/regr_pga_Mc8.5_0.0_6.5_VR_99.4.pckl'
#
## REad in the fixed model to get coefficients:
#invfile = open(fixedmodel,'r')
#inv_dat = pickle.load(invfile)
#invfile.close()

print 'Starting iterative mixed effects...'
# Fix which coefficients?
a2 = inv_dat.m[1]
a3 = inv_dat.m[2]
a4 = inv_dat.m[3]
a5 = inv_dat.m[4]
a6 = 'none'

vs30_correct = 0

fixedcoeff_list = 'a2_%.2f_a3_%.2f_a4_%.2f_a5_%.2f' % (a2,a3,a4,a5)

print 'Running mixed effects'
#Now try with mixed effects:
run_name = 'mixedregr_v3anza2013_pga_%icoeff_%s_Mc_8.5_res4_noVs30' % (ncoeff,fixedcoeff_list)
run_home=home+'/models/residuals/'
resaxlim_r = [[0,180],[-5,5]]
resaxlim_mw = [[0,4],[-5,5]]



#Fictitious depth parameter:
c=4.5

#Initialize the residuals directories:
inithome=HOME+'/anza/models/residuals/'

runall=1

#Initialize directories:
runall=run_res.init(inithome,run_name)

if runall==0:
    print 'Not clobbering, exiting...'
    
elif runall==1:
    print 'Continuing...'

    # Now run mixed effects approach #
    invdat,invpath,tresid,mixed_residuals,d_r_prediction,mixed_resid_path=run_inv.run_mixedeffects(home,codehome,run_name,dbpath,dbname,Mc,vref,c,predictive_parameter='pga',ncoeff=ncoeff,data_correct=vs30_correct,a1='none',a2=a2,a3=a3,a4=a4,a5=a5,a6=a6)
    
    
    print 'Variance reduction is %f' % invdat.VR
    
    print 'Plotting mixed effects model with data'
    # Plot data with model:
    mixedinv = run_inv.plot_data_model(home,dbpath,dbname,invpath,coeff_file,mdep_ffdf,plotdist,ask_dist,Mc,axlims,bmin,bmax,vref,predictive_parameter=predictive_parameter,ncoeff=ncoeff,data_correct=vs30_correct)

    ## If just plotting:
    #invpath=setmixedmodel
    #mixedinv = run_inv.plot_data_model(home,dbpath,dbname,invpath,coeff_file,mdep_ffdf,plotdist,ask_dist,Mc,axlims,bmin,bmax,vref)


    # Get some mean values for the stats file:
    mean_tot=np.mean(mixed_residuals.total_residual)
    std_dev_tot=np.std(mixed_residuals.total_residual)
    
    E_mean=mixed_residuals.E_mean
    E_std_dev=mixed_residuals.E_std

    print 'plotting all residuals'
    # Plot all residuals:
    run_res.plot_total(tresid,home,run_name,resaxlim_mw)

    print 'plotting event terms'
    # Plot event terms:
    eventfig1 = run_res.makeEvents_mixed(run_home,run_name,mixed_resid_path,Mc,vref,mdep_ffdf,resaxlim_mw)

    # Get station objects:
    station_list=run_res.makeStations_mixed(run_home,run_name,mixed_resid_path)

    # Plot W residuals on one plot:
    W_mean,W_std_dev = run_res.plot_Wresid(run_home,run_name,resaxlim_mw)

    print 'plotting station terms'
    # Plot by station:
    run_res.plot_site_WE(run_home,run_name,resaxlim_mw)

    print 'plotting path residuals'
    # Plot path residuals:
    f_mw,f_dist,pterm_mean,pterm_std=run_res.plotPath_mixed(run_home,run_name,mixed_resid_path,resaxlim_mw,resaxlim_r)
    
    print 'Writing mixed stats to a file'
    # Writing stats to a file
    run_res.write_stats_mixed(run_home,run_name,mixed_resid_path,mean_tot,std_dev_tot,E_mean,E_std_dev,W_mean,W_std_dev,pterm_mean,pterm_std)
    
    

###################################################################################
#################       MIXED EFFECTS with a4 and a5 fixed      ###################
###################################################################################
#

# Read in model so you don't have to run the entire thing:
fixedmodel = '/Users/vsahakian/anza/models/pckl/v3anza2013_pga_noVs30_5coeff/regr_pga_Mc8.5_0.0_6.5_VR_99.3.pckl'

# REad in the fixed model to get coefficients:
invfile = open(fixedmodel,'r')
inv_dat = pickle.load(invfile)
invfile.close()

print 'Starting iterative mixed effects...'
# Fix which coefficients?
a2 = 'none'
a3 = 'none'
a4 = inv_dat.m[3]
a5 = inv_dat.m[4]
a6 = 'none'

vs30_correct = 0

fixedcoeff_list = 'a4_%.2f_a5_%.2f' % (a4,a5)

print 'Running mixed effects'
#Now try with mixed effects:
run_name = 'mixedregr_v3anza2013_pga_%icoeff_%s_Mc_8.5_res4_noVs30' % (ncoeff,fixedcoeff_list)
run_home=home+'/models/residuals/'
resaxlim_r = [[0,180],[-5,5]]
resaxlim_mw = [[0,4],[-5,5]]

# INversion run  name - called "dbname" here - setting a4 through a6:
dbname = 'v3anza2013_pga_noVs30_5coeff_a4_a5'

#Fictitious depth parameter:
c=4.5

#Initialize the residuals directories:
inithome=HOME+'/anza/models/residuals/'

runall=1

#Initialize directories:
runall=run_res.init(inithome,run_name)

if runall==0:
    print 'Not clobbering, exiting...'
    
elif runall==1:
    print 'Continuing...'

    # Now run mixed effects approach #
    invdat,invpath,tresid,mixed_residuals,d_r_prediction,mixed_resid_path=run_inv.run_mixedeffects(home,codehome,run_name,dbpath,dbname,Mc,vref,c,predictive_parameter='pga',ncoeff=ncoeff,data_correct=vs30_correct,a1='none',a2=a2,a3=a3,a4=a4,a5=a5,a6=a6)
    
    
    print 'Variance reduction is %f' % invdat.VR
    
    print 'Plotting mixed effects model with data'
    # Plot data with model:
    mixedinv = run_inv.plot_data_model(home,dbpath,dbname,invpath,coeff_file,mdep_ffdf,plotdist,ask_dist,Mc,axlims,bmin,bmax,vref,predictive_parameter=predictive_parameter,ncoeff=ncoeff,data_correct=vs30_correct)

    ## If just plotting:
    #invpath=setmixedmodel
    #mixedinv = run_inv.plot_data_model(home,dbpath,dbname,invpath,coeff_file,mdep_ffdf,plotdist,ask_dist,Mc,axlims,bmin,bmax,vref)


    # Get some mean values for the stats file:
    mean_tot=np.mean(mixed_residuals.total_residual)
    std_dev_tot=np.std(mixed_residuals.total_residual)
    
    E_mean=mixed_residuals.E_mean
    E_std_dev=mixed_residuals.E_std

    print 'plotting all residuals'
    # Plot all residuals:
    run_res.plot_total(tresid,home,run_name,resaxlim_mw)

    print 'plotting event terms'
    # Plot event terms:
    eventfig1 = run_res.makeEvents_mixed(run_home,run_name,mixed_resid_path,Mc,vref,mdep_ffdf,resaxlim_mw)

    # Get station objects:
    station_list=run_res.makeStations_mixed(run_home,run_name,mixed_resid_path)

    # Plot W residuals on one plot:
    W_mean,W_std_dev = run_res.plot_Wresid(run_home,run_name,resaxlim_mw)

    print 'plotting station terms'
    # Plot by station:
    run_res.plot_site_WE(run_home,run_name,resaxlim_mw)

    print 'plotting path residuals'
    # Plot path residuals:
    f_mw,f_dist,pterm_mean,pterm_std=run_res.plotPath_mixed(run_home,run_name,mixed_resid_path,resaxlim_mw,resaxlim_r)
    
    print 'Writing mixed stats to a file'
    # Writing stats to a file
    run_res.write_stats_mixed(run_home,run_name,mixed_resid_path,mean_tot,std_dev_tot,E_mean,E_std_dev,W_mean,W_std_dev,pterm_mean,pterm_std)

    
        
                
###################################################################################
#################       MIXED EFFECTS with a4 fixed to 1.0      ###################
###################################################################################
#

# Read in model so you don't have to run the entire thing:
fixedmodel = '/Users/vsahakian/anza/models/pckl/v3anza2013_pga_noVs30_5coeff/regr_pga_Mc8.5_0.0_6.5_VR_99.3.pckl'

# REad in the fixed model to get coefficients:
invfile = open(fixedmodel,'r')
inv_dat = pickle.load(invfile)
invfile.close()

print 'Starting iterative mixed effects...'
# Fix which coefficients?
a2 = 'none'
a3 = 'none'
a4 = -1.0
a5 = 'none'
a6 = 'none'

vs30_correct = 0

fixedcoeff_list = 'a4_%.2f' % (a4)

print 'Running mixed effects'
#Now try with mixed effects:
run_name = 'mixedregr_v3anza2013_pga_%icoeff_%s_Mc_8.5_res4_noVs30' % (ncoeff,fixedcoeff_list)
run_home=home+'/models/residuals/'
resaxlim_r = [[0,180],[-5,5]]
resaxlim_mw = [[0,4],[-5,5]]

# INversion run  name - called "dbname" here - setting a4 through a6:
dbname = 'v3anza2013_pga_noVs30_5coeff_a4_-1.0'

#Fictitious depth parameter:
c=4.5

#Initialize the residuals directories:
inithome=HOME+'/anza/models/residuals/'

runall=1

#Initialize directories:
runall=run_res.init(inithome,run_name)

if runall==0:
    print 'Not clobbering, exiting...'
    
elif runall==1:
    print 'Continuing...'

    # Now run mixed effects approach #
    invdat,invpath,tresid,mixed_residuals,d_r_prediction,mixed_resid_path=run_inv.run_mixedeffects(home,codehome,run_name,dbpath,dbname,Mc,vref,c,predictive_parameter='pga',ncoeff=ncoeff,data_correct=vs30_correct,a1='none',a2=a2,a3=a3,a4=a4,a5=a5,a6=a6)
    
    
    print 'Variance reduction is %f' % invdat.VR
    
    print 'Plotting mixed effects model with data'
    # Plot data with model:
    mixedinv = run_inv.plot_data_model(home,dbpath,dbname,invpath,coeff_file,mdep_ffdf,plotdist,ask_dist,Mc,axlims,bmin,bmax,vref,predictive_parameter=predictive_parameter,ncoeff=ncoeff,data_correct=vs30_correct)

    ## If just plotting:
    #invpath=setmixedmodel
    #mixedinv = run_inv.plot_data_model(home,dbpath,dbname,invpath,coeff_file,mdep_ffdf,plotdist,ask_dist,Mc,axlims,bmin,bmax,vref)


    # Get some mean values for the stats file:
    mean_tot=np.mean(mixed_residuals.total_residual)
    std_dev_tot=np.std(mixed_residuals.total_residual)
    
    E_mean=mixed_residuals.E_mean
    E_std_dev=mixed_residuals.E_std

    print 'plotting all residuals'
    # Plot all residuals:
    run_res.plot_total(tresid,home,run_name,resaxlim_mw)

    print 'plotting event terms'
    # Plot event terms:
    eventfig1 = run_res.makeEvents_mixed(run_home,run_name,mixed_resid_path,Mc,vref,mdep_ffdf,resaxlim_mw)

    # Get station objects:
    station_list=run_res.makeStations_mixed(run_home,run_name,mixed_resid_path)

    # Plot W residuals on one plot:
    W_mean,W_std_dev = run_res.plot_Wresid(run_home,run_name,resaxlim_mw)

    print 'plotting station terms'
    # Plot by station:
    run_res.plot_site_WE(run_home,run_name,resaxlim_mw)

    print 'plotting path residuals'
    # Plot path residuals:
    f_mw,f_dist,pterm_mean,pterm_std=run_res.plotPath_mixed(run_home,run_name,mixed_resid_path,resaxlim_mw,resaxlim_r)
    
    print 'Writing mixed stats to a file'
    # Writing stats to a file
    run_res.write_stats_mixed(run_home,run_name,mixed_resid_path,mean_tot,std_dev_tot,E_mean,E_std_dev,W_mean,W_std_dev,pterm_mean,pterm_std)
    



###################################################################################
#################       MIXED EFFECTS with a4 fixed to 1.2      ###################
###################################################################################
#

# Read in model so you don't have to run the entire thing:
fixedmodel = '/Users/vsahakian/anza/models/pckl/v3anza2013_pga_noVs30_5coeff/regr_pga_Mc8.5_0.0_6.5_VR_99.3.pckl'

# REad in the fixed model to get coefficients:
invfile = open(fixedmodel,'r')
inv_dat = pickle.load(invfile)
invfile.close()

print 'Starting iterative mixed effects...'
# Fix which coefficients?
a2 = 'none'
a3 = 'none'
a4 = -1.2
a5 = 'none'
a6 = 'none'

vs30_correct = 0

fixedcoeff_list = 'a4_%.2f' % (a4)

print 'Running mixed effects'
#Now try with mixed effects:
run_name = 'mixedregr_v3anza2013_pga_%icoeff_%s_Mc_8.5_res4_noVs30' % (ncoeff,fixedcoeff_list)
run_home=home+'/models/residuals/'
resaxlim_r = [[0,180],[-5,5]]
resaxlim_mw = [[0,4],[-5,5]]

# INversion run  name - called "dbname" here - setting a4 through a6:
dbname = 'v3anza2013_pga_noVs30_5coeff_a4_-1.2'

#Fictitious depth parameter:
c=4.5

#Initialize the residuals directories:
inithome=HOME+'/anza/models/residuals/'

runall=1

#Initialize directories:
runall=run_res.init(inithome,run_name)

if runall==0:
    print 'Not clobbering, exiting...'
    
elif runall==1:
    print 'Continuing...'

    # Now run mixed effects approach #
    invdat,invpath,tresid,mixed_residuals,d_r_prediction,mixed_resid_path=run_inv.run_mixedeffects(home,codehome,run_name,dbpath,dbname,Mc,vref,c,predictive_parameter='pga',ncoeff=ncoeff,data_correct=vs30_correct,a1='none',a2=a2,a3=a3,a4=a4,a5=a5,a6=a6)
    
    
    print 'Variance reduction is %f' % invdat.VR
    
    print 'Plotting mixed effects model with data'
    # Plot data with model:
    mixedinv = run_inv.plot_data_model(home,dbpath,dbname,invpath,coeff_file,mdep_ffdf,plotdist,ask_dist,Mc,axlims,bmin,bmax,vref,predictive_parameter=predictive_parameter,ncoeff=ncoeff,data_correct=vs30_correct)

    ## If just plotting:
    #invpath=setmixedmodel
    #mixedinv = run_inv.plot_data_model(home,dbpath,dbname,invpath,coeff_file,mdep_ffdf,plotdist,ask_dist,Mc,axlims,bmin,bmax,vref)


    # Get some mean values for the stats file:
    mean_tot=np.mean(mixed_residuals.total_residual)
    std_dev_tot=np.std(mixed_residuals.total_residual)
    
    E_mean=mixed_residuals.E_mean
    E_std_dev=mixed_residuals.E_std

    print 'plotting all residuals'
    # Plot all residuals:
    run_res.plot_total(tresid,home,run_name,resaxlim_mw)

    print 'plotting event terms'
    # Plot event terms:
    eventfig1 = run_res.makeEvents_mixed(run_home,run_name,mixed_resid_path,Mc,vref,mdep_ffdf,resaxlim_mw)

    # Get station objects:
    station_list=run_res.makeStations_mixed(run_home,run_name,mixed_resid_path)

    # Plot W residuals on one plot:
    W_mean,W_std_dev = run_res.plot_Wresid(run_home,run_name,resaxlim_mw)

    print 'plotting station terms'
    # Plot by station:
    run_res.plot_site_WE(run_home,run_name,resaxlim_mw)

    print 'plotting path residuals'
    # Plot path residuals:
    f_mw,f_dist,pterm_mean,pterm_std=run_res.plotPath_mixed(run_home,run_name,mixed_resid_path,resaxlim_mw,resaxlim_r)
    
    print 'Writing mixed stats to a file'
    # Writing stats to a file
    run_res.write_stats_mixed(run_home,run_name,mixed_resid_path,mean_tot,std_dev_tot,E_mean,E_std_dev,W_mean,W_std_dev,pterm_mean,pterm_std)



###################################################################################
#################       MIXED EFFECTS with a4 fixed to 1.5      ###################
###################################################################################
#

# Read in model so you don't have to run the entire thing:
fixedmodel = '/Users/vsahakian/anza/models/pckl/v3anza2013_pga_noVs30_5coeff/regr_pga_Mc8.5_0.0_6.5_VR_99.3.pckl'

# REad in the fixed model to get coefficients:
invfile = open(fixedmodel,'r')
inv_dat = pickle.load(invfile)
invfile.close()

print 'Starting iterative mixed effects...'
# Fix which coefficients?
a2 = 'none'
a3 = 'none'
a4 = -1.5
a5 = 'none'
a6 = 'none'

vs30_correct = 0

fixedcoeff_list = 'a4_%.2f' % (a4)

print 'Running mixed effects'
#Now try with mixed effects:
run_name = 'mixedregr_v3anza2013_pga_%icoeff_%s_Mc_8.5_res4_noVs30' % (ncoeff,fixedcoeff_list)
run_home=home+'/models/residuals/'
resaxlim_r = [[0,180],[-5,5]]
resaxlim_mw = [[0,4],[-5,5]]

# INversion run  name - called "dbname" here - setting a4 through a6:
dbname = 'v3anza2013_pga_noVs30_5coeff_a4_-1.5'

#Fictitious depth parameter:
c=4.5

#Initialize the residuals directories:
inithome=HOME+'/anza/models/residuals/'

runall=1

#Initialize directories:
runall=run_res.init(inithome,run_name)

if runall==0:
    print 'Not clobbering, exiting...'
    
elif runall==1:
    print 'Continuing...'

    # Now run mixed effects approach #
    invdat,invpath,tresid,mixed_residuals,d_r_prediction,mixed_resid_path=run_inv.run_mixedeffects(home,codehome,run_name,dbpath,dbname,Mc,vref,c,predictive_parameter='pga',ncoeff=ncoeff,data_correct=vs30_correct,a1='none',a2=a2,a3=a3,a4=a4,a5=a5,a6=a6)
    
    
    print 'Variance reduction is %f' % invdat.VR
    
    print 'Plotting mixed effects model with data'
    # Plot data with model:
    mixedinv = run_inv.plot_data_model(home,dbpath,dbname,invpath,coeff_file,mdep_ffdf,plotdist,ask_dist,Mc,axlims,bmin,bmax,vref,predictive_parameter=predictive_parameter,ncoeff=ncoeff,data_correct=vs30_correct)

    ## If just plotting:
    #invpath=setmixedmodel
    #mixedinv = run_inv.plot_data_model(home,dbpath,dbname,invpath,coeff_file,mdep_ffdf,plotdist,ask_dist,Mc,axlims,bmin,bmax,vref)


    # Get some mean values for the stats file:
    mean_tot=np.mean(mixed_residuals.total_residual)
    std_dev_tot=np.std(mixed_residuals.total_residual)
    
    E_mean=mixed_residuals.E_mean
    E_std_dev=mixed_residuals.E_std

    print 'plotting all residuals'
    # Plot all residuals:
    run_res.plot_total(tresid,home,run_name,resaxlim_mw)

    print 'plotting event terms'
    # Plot event terms:
    eventfig1 = run_res.makeEvents_mixed(run_home,run_name,mixed_resid_path,Mc,vref,mdep_ffdf,resaxlim_mw)

    # Get station objects:
    station_list=run_res.makeStations_mixed(run_home,run_name,mixed_resid_path)

    # Plot W residuals on one plot:
    W_mean,W_std_dev = run_res.plot_Wresid(run_home,run_name,resaxlim_mw)

    print 'plotting station terms'
    # Plot by station:
    run_res.plot_site_WE(run_home,run_name,resaxlim_mw)

    print 'plotting path residuals'
    # Plot path residuals:
    f_mw,f_dist,pterm_mean,pterm_std=run_res.plotPath_mixed(run_home,run_name,mixed_resid_path,resaxlim_mw,resaxlim_r)
    
    print 'Writing mixed stats to a file'
    # Writing stats to a file
    run_res.write_stats_mixed(run_home,run_name,mixed_resid_path,mean_tot,std_dev_tot,E_mean,E_std_dev,W_mean,W_std_dev,pterm_mean,pterm_std)
    

###################################################################################
#################       MIXED EFFECTS with a5 fixed to 0.02     ###################
###################################################################################
#

# Read in model so you don't have to run the entire thing:
fixedmodel = '/Users/vsahakian/anza/models/pckl/v3anza2013_pga_noVs30_5coeff/regr_pga_Mc8.5_0.0_6.5_VR_99.3.pckl'

# REad in the fixed model to get coefficients:
invfile = open(fixedmodel,'r')
inv_dat = pickle.load(invfile)
invfile.close()

print 'Starting iterative mixed effects...'
# Fix which coefficients?
a2 = 'none'
a3 = 'none'
a4 = 'none'
a5 = -0.02
a6 = 'none'

vs30_correct = 0

fixedcoeff_list = 'a5_%.2f' % (a5)

print 'Running mixed effects'
#Now try with mixed effects:
run_name = 'mixedregr_v3anza2013_pga_%icoeff_%s_Mc_8.5_res4_noVs30' % (ncoeff,fixedcoeff_list)
run_home=home+'/models/residuals/'
resaxlim_r = [[0,180],[-5,5]]
resaxlim_mw = [[0,4],[-5,5]]

# INversion run  name - called "dbname" here - setting a4 through a6:
dbname = 'v3anza2013_pga_noVs30_5coeff_a5_-0.02'

#Fictitious depth parameter:
c=4.5

#Initialize the residuals directories:
inithome=HOME+'/anza/models/residuals/'

runall=1

#Initialize directories:
runall=run_res.init(inithome,run_name)

if runall==0:
    print 'Not clobbering, exiting...'
    
elif runall==1:
    print 'Continuing...'

    # Now run mixed effects approach #
    invdat,invpath,tresid,mixed_residuals,d_r_prediction,mixed_resid_path=run_inv.run_mixedeffects(home,codehome,run_name,dbpath,dbname,Mc,vref,c,predictive_parameter='pga',ncoeff=ncoeff,data_correct=vs30_correct,a1='none',a2=a2,a3=a3,a4=a4,a5=a5,a6=a6)
    
    
    print 'Variance reduction is %f' % invdat.VR
    
    print 'Plotting mixed effects model with data'
    # Plot data with model:
    mixedinv = run_inv.plot_data_model(home,dbpath,dbname,invpath,coeff_file,mdep_ffdf,plotdist,ask_dist,Mc,axlims,bmin,bmax,vref,predictive_parameter=predictive_parameter,ncoeff=ncoeff,data_correct=vs30_correct)

    ## If just plotting:
    #invpath=setmixedmodel
    #mixedinv = run_inv.plot_data_model(home,dbpath,dbname,invpath,coeff_file,mdep_ffdf,plotdist,ask_dist,Mc,axlims,bmin,bmax,vref)


    # Get some mean values for the stats file:
    mean_tot=np.mean(mixed_residuals.total_residual)
    std_dev_tot=np.std(mixed_residuals.total_residual)
    
    E_mean=mixed_residuals.E_mean
    E_std_dev=mixed_residuals.E_std

    print 'plotting all residuals'
    # Plot all residuals:
    run_res.plot_total(tresid,home,run_name,resaxlim_mw)

    print 'plotting event terms'
    # Plot event terms:
    eventfig1 = run_res.makeEvents_mixed(run_home,run_name,mixed_resid_path,Mc,vref,mdep_ffdf,resaxlim_mw)

    # Get station objects:
    station_list=run_res.makeStations_mixed(run_home,run_name,mixed_resid_path)

    # Plot W residuals on one plot:
    W_mean,W_std_dev = run_res.plot_Wresid(run_home,run_name,resaxlim_mw)

    print 'plotting station terms'
    # Plot by station:
    run_res.plot_site_WE(run_home,run_name,resaxlim_mw)

    print 'plotting path residuals'
    # Plot path residuals:
    f_mw,f_dist,pterm_mean,pterm_std=run_res.plotPath_mixed(run_home,run_name,mixed_resid_path,resaxlim_mw,resaxlim_r)
    
    print 'Writing mixed stats to a file'
    # Writing stats to a file
    run_res.write_stats_mixed(run_home,run_name,mixed_resid_path,mean_tot,std_dev_tot,E_mean,E_std_dev,W_mean,W_std_dev,pterm_mean,pterm_std)