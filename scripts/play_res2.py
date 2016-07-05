######Residual Computatoin######

import run_res

home='/Users/vsahakian/anza/models/residuals/'
run_name='abdb_0-6.5'
dbpath='/Users/vsahakian/anza/data/abdb.pckl'
modelpath='/Users/vsahakian/anza/models/pckl/regr_0.0_6.5_resid_2676.06963031.pckl'

ffdf_flag=0
resaxlim=[[1,4],[-4,4]]

#Initialize directories:
run_res.init(home,run_name)

#Get total residuals and plots:
tr_mw,tot_resid=run_res.get_total_res(home,run_name,dbpath,modelpath,ffdf_flag,resaxlim)

#Get event/within-event residuals:
E_evnum,E_mw,E_residual,E_mean,E_std_dev=run_res.getEW_makeEvents(home,run_name,dbpath,modelpath,ffdf_flag)

#Make station objects:
run_res.sta_list(home,run_name,dbpath)

#Plot within-event residuals by station, and save to file:
run_res.plot_Wresid(home,run_name,resaxlim)