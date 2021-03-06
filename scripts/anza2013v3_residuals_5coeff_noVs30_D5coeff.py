######Residual Computatoin######

import run_res
from os import path
import raytracing as rt

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
    


home=HOME+'/anza/models/residuals/'
run_name='v3anza2013_Mc8.5_pgrid_5sta_res4_pga_noVs30_5coeff_D5coeff'
dbpath=HOME+'/anza/data/databases/v3anza2013/v3anza2013_pgrid_5sta_res4.pckl'
modelpath=HOME+'/anza/models/pckl/v3anza2013_pga_noVs30_5coeff_a4_a5/mixedregr_v3anza2013_pga_noVs30_5coeff_a4_a5_pga__ncoeff5_Mc_8.5_VR_99.4_a4_-1.72_a5_-0.01.pckl'
#rayfile_vp='/media/vsahakian/katmai/anza/fm3d/abdb_5sta_topography/Vp/rays.dat'
#rayfile_vs='/media/vsahakian/katmai/anza/fm3d/abdb_5sta_topography/Vs/rays.dat'

vref=760


########

pred_param='pga'
ncoeff=5
vs30_correct=0

ffdf_flag=0
Mc=8.5
resaxlim=[[0,5],[-4.5,4.5]]
resaxlim_dist=[[0,190],[-4.5,4.5]]

#Initialize runall as being 1, in case the directory is new:
runall=1

#Initialize directories:
runall=run_res.init(home,run_name)

if runall==0:
    print 'Not clobbering, exiting...'
    
elif runall==1:
    print 'Continuing...'
    
    #Get total residuals and plots:
    tr_mw,tot_resid,mean_tot,std_dev_tot=run_res.get_total_res(home,run_name,dbpath,modelpath,Mc,ffdf_flag,resaxlim,predictive_parameter=pred_param,ncoeff=ncoeff,data_correct=vs30_correct)
    
    #Get event/within-event residuals:
    E_evnum,E_mw,E_residual,E_mean,E_std_dev=run_res.getEW_makeEvents(home,run_name,dbpath,modelpath,Mc,vref,ffdf_flag,resaxlim,predictive_parameter=pred_param,ncoeff=ncoeff,data_correct=vs30_correct)
    
    #Make station objects:
    run_res.sta_list(home,run_name,dbpath)
    
    #Plot within-event residuals by station on one plot, and save to file:
    W_mean,W_std_dev=run_res.plot_Wresid(home,run_name,resaxlim)
    
    #Plot within-event residuals by station in subplots:
    run_res.plot_site_WE(home,run_name,resaxlim)
    
    #Write to an overall residuals plus more object, get path residual...
    f_mw,f_dist,allresiduals,pterm_mean,pterm_std=run_res.get_path_resid_make_object(home,run_name,dbpath,resaxlim,resaxlim_dist)
    
    #Put this stuff in play res
    ##Read in the Vp and Vs ray info from teh rayfile:
    #vp_path_list,rec_id,src_id=rt.parse_rayfile(rayfile_vp)
    #vs_path_list,rec_id,src_id=rt.parse_rayfile(rayfile_vs)
    #
    ##Store the ray info in the residuals object:
    #allresiduals
    
    #Write to file:
    run_res.write_stats(home,run_name,mean_tot,std_dev_tot,E_mean,E_std_dev,W_mean,W_std_dev,pterm_mean,pterm_std)