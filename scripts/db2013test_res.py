######Residual Computatoin######

import run_res
from os import path
import raytracing as rt

#Change this parameter depending on where you run:
#0=desktop
#1=mac

what_home=0

if what_home==0:
    #Desktop:
    HOME='/media/vsahakian/katmai'
elif what_home==1:
    #Mac:
    HOME='/Users/vsahakian'
    

##Centered at M 8.5 
#home=HOME+'/anza/models/residuals/'
#run_name='db2013test_0-6.5'
#dbpath=HOME+'/anza/data/databases/db2013_test/db2013test_5sta.pckl'
#modelpath=HOME+'/anza/models/pckl/test2013/regr_0.0_6.5_VR_99.2.pckl'

##Centered at M 3.0
#home=HOME+'/anza/models/residuals/'
#run_name='db2013test_0-6.5_Mc3'
#dbpath=HOME+'/anza/data/databases/db2013_test/db2013test_5sta.pckl'
#modelpath=HOME+'/anza/models/pckl/test2013/regr_Mc3_0.0_6.5_VR_99.2.pckl'

##Centered at M 7.0
#home=HOME+'/anza/models/residuals/'
#run_name='db2013test_0-6.5_Mc7'
#dbpath=HOME+'/anza/data/databases/db2013_test/db2013test_5sta.pckl'
#modelpath=HOME+'/anza/models/pckl/test2013/regr_Mc7_0.0_6.5_VR_99.2.pckl'

##Centered at M 9.0
#home=HOME+'/anza/models/residuals/'
#run_name='db2013test_0-6.5_Mc9'
#dbpath=HOME+'/anza/data/databases/db2013_test/db2013test_5sta.pckl'
#modelpath=HOME+'/anza/models/pckl/test2013/regr_Mc9_0.0_6.5_VR_99.2.pckl'

##Centered at M 8.0
#home=HOME+'/anza/models/residuals/'
#run_name='db2013test_0-6.5_Mc8'
#dbpath=HOME+'/anza/data/databases/db2013_test/db2013test_5sta.pckl'
#modelpath=HOME+'/anza/models/pckl/test2013/regr_Mc8_0.0_6.5_VR_99.2.pckl'

#Centered at M 8.1
home=HOME+'/anza/models/residuals/'
run_name='db2013test_0-6.5_Mc8.1'
dbpath=HOME+'/anza/data/databases/db2013_test/db2013test_5sta.pckl'
modelpath=HOME+'/anza/models/pckl/test2013/regr_Mc8.1_0.0_6.5_VR_99.2.pckl'

#rayfile_vp='/media/vsahakian/katmai/anza/fm3d/abdb_5sta_topography/Vp/rays.dat'
#rayfile_vs='/media/vsahakian/katmai/anza/fm3d/abdb_5sta_topography/Vs/rays.dat'


########

ffdf_flag=0
resaxlim=[[0,4],[-4,4]]
resaxlim_dist=[[2,260],[-4,4]]

#Initialize runall as being 1, in case the directory is new:
runall=1

#Initialize directories:
runall=run_res.init(home,run_name)

if runall==0:
    print 'Not clobbering, exiting...'
    
elif runall==1:
    print 'Continuing...'
    
    #Get total residuals and plots:
    tr_mw,tot_resid,mean_tot,std_dev_tot=run_res.get_total_res(home,run_name,dbpath,modelpath,ffdf_flag,resaxlim)
    
    #Get event/within-event residuals:
    E_evnum,E_mw,E_residual,E_mean,E_std_dev=run_res.getEW_makeEvents(home,run_name,dbpath,modelpath,ffdf_flag,resaxlim)
    
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