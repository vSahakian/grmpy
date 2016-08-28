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
    

#home=HOME+'/anza/models/residuals/'
#run_name='abdb_0-6.5_addindex'
#dbpath=HOME+'/anza/data/abdb.pckl'
#modelpath=HOME+'/anza/models/pckl/regr_0.0_6.5_resid_2676.06963031.pckl'
#
#home=HOME+'/anza/models/residuals/'
#run_name='abdb_0-2-3-4-6.5'
#dbpath=HOME+'/anza/data/abdb.pckl'
#modelpath=HOME+'/anza/models/pckl/regr_0.0_2.0_3.0_4.0_6.5_VR_98.9.pckl'

#home=HOME+'/anza/models/residuals/'
#run_name='abdb_0-3-6.5'
#dbpath=HOME+'/anza/data/abdb.pckl'
#modelpath=HOME+'/anza/models/pckl/regr_0.0_3.0_6.5_VR_98.9.pckl'

#home=HOME+'/anza/models/residuals/'
#run_name='abdb_0-1-2-6.5'
#dbpath=HOME+'/anza/data/abdb.pckl'
#modelpath=HOME+'/anza/models/pckl/regr_0.0_1.0_2.0_6.5_VR_98.9.pckl'

#home=HOME+'/anza/models/residuals/'
#run_name='abdb_0-6.5_VR'
#dbpath=HOME+'/anza/data/abdb.pckl'
#modelpath=HOME+'/anza/models/pckl/regr_0.0_6.5_VR_98.9.pckl'

home=HOME+'/anza/models/residuals/'
run_name='abdb_5sta_0-6.5_VR'
dbpath=HOME+'/anza/data/abdb_5sta.pckl'
modelpath=HOME+'/anza/models/pckl/regr_0.0_6.5_VR_98.9.pckl'
rayfile_vp='/Users/vsahakian/anza/data/vm/fulltest_Vp/rays.dat'
rayfile_vs='/Users/vsahakian/anza/data/vm/fulltest_Vs/rays.dat'



########

ffdf_flag=0
resaxlim=[[1,4],[-4,4]]

#Initialize runall as being 1, in case the directory is new:
runall=1

#Initialize directories:
runall=run_res.init(home,run_name)

if runall==0:
    print 'Not clobbering, exiting...'
    
elif runall==1:
    print 
    
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
    allresiduals=run_res.get_path_resid_make_object(home,run_name,dbpath)
    
    #Put this stuff in play res
    ##Read in the Vp and Vs ray info from teh rayfile:
    #vp_path_list,rec_id,src_id=rt.parse_rayfile(rayfile_vp)
    #vs_path_list,rec_id,src_id=rt.parse_rayfile(rayfile_vs)
    #
    ##Store the ray info in the residuals object:
    #allresiduals
    
    #Write to file:
    run_res.write_stats(home,run_name,mean_tot,std_dev_tot,E_mean,E_std_dev,W_mean,W_std_dev)