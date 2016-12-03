##Play around iwth the raytracing...
#VJS 8/2016

import run_res
from os import path
import dread
import cPickle as pickle
import raytracing as rt
import matplotlib.pyplot as plt

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

home=HOME+'/anza/models/residuals/'
#run_name='anza2013_Mc8.5_pgrid_5sta'
run_name='mixedregr_anza2013_Mc_8.5'
dbpath=HOME+'/anza/data/databases/anza2013/anza2013_pgrid_5sta.pckl'
modelpath=HOME+'/anza/models/pckl/anza2013/regr_Mc8.5_0.0_6.5_VR_99.4.pckl'
residname='mixedregr_anza2013_Mc_8.5_VR_99.9'
#rayfile_vp='/media/vsahakian/katmai/anza/fm3d/anza2013_Mc8.5_pgrid_5sta/Vp/rays.dat'
rayfile_vs='/media/vsahakian/katmai/anza/fm3d/anza2013_Mc8.5_pgrid_5sta/Vs/rays.dat'
faultfile='/media/vsahakian/katmai/anza/data/faults/Holocene_LatestPleistocene_117.5w_115.5w_33n_34n.pckl'

#########################
#####Set up Raytracing###
#########################
#
##Set up the sources.in and receivers.in files for Vs.  So far, Vp has been run:
##in /media/vsahakian/katmai/anza/data/vm/fulltest_Vp
#
##Run these with raytracing write_sourcein and write_receiverin
##Use the directories defined in lines 15 - 33
##Import raytracing above
##Set velocity type:
#veltype=1
##Set longitude format (here, 241 instead of -118 for lontype=1)
#lontype=1  #for Malcolm's new propgrid format
#
##Write source.in file:
#rt.write_sourcein(home,run_name,veltype,lontype)
#
##Write recievers.in file:
#rt.write_receiverin(home,run_name,lontype)




#####################
###Store ray files###
#####################

##Get the "in" and "out" residual object file names:
##"In" is the original residuals object:
#
#Get the run directory:
run_dir=path.expanduser(home+run_name+'/')
#Get the residuals object:
#residfile_in=run_dir+run_name+'_robj.pckl'
residfile_in=run_dir+run_name+'_VR_99.9_robj.pckl'

#"Out" is the _raydat.pckl object - also serves as the "in" for the second:
rbase=residfile_in.split('.pckl')
#residfile_out=rbase[0]+'_raydat.pckl'
residfile_out=rbase[0]+'_raydat.pckl'

##Longitude conversion?  
#Yes - these rayfiles are in positive west, convert to negative west (so flag =1)
lonconvert=1

##For vp:
#rayfile=rayfile_vp
#veltype=1
#
##Read in:
#rt.store_rayinfo(residfile_in,residfile_out,rayfile,veltype,lonconvert)
#print 'Stored Vp rays in the residuals object'

#####

#For vs:
rayfile=rayfile_vs
veltype=2

##Read in:
#rt.store_rayinfo(residfile_out,residfile_out,rayfile,veltype,lonconvert)
#print 'Stored Vs rays in the residuals object'

#Read in:
rt.store_rayinfo(residfile_in,residfile_out,rayfile,veltype,lonconvert)
print 'Stored Vs rays in the residuals object'




###############
###Plot rays###
###############


#Vp and Vs:
veltype=[1,2]

#map, and cross sections:
view=[0,1,2]
axlims_view0=[[-118.0,-115.2],[32.3,34.5]]
axlims_view1=[[32.3,34.5],[-28,5]]
axlims_view2=[[-118.0,-115.2],[-28,5]]
axlims=[axlims_view0,axlims_view1,axlims_view2]

stations=1
events=1
by_path=1
mymap='jet'
cutoff_val=1.0


###
#Plot in a loop
#for vel_i in range(len(veltype)):
#    for view_i in range(len(view)):
#        rt.plot_rays(home,run_name,residname,veltype[vel_i],view[view_i],axlims[view_i],stations,events,by_path,mymap,faultfile)
#        print 'Plotted velocity type '+str(veltype[vel_i])+', view '+str(view[view_i])
#        plt.close()
#        plt.close()

##Plot with the cutoff value:
##Plot in a loop
#for vel_i in range(len(veltype)):
#    for view_i in range(len(view)):
#        rt.plot_rays_cutoffval(home,run_name,residname,veltype[vel_i],view[view_i],axlims[view_i],stations,events,mymap,faultfile,cutoff_val)
#        print 'Plotted cuttoff value '+str(cutoff_val)+', velocity type '+str(veltype[vel_i])+', view '+str(view[view_i])
#        plt.close()
#        plt.close()
        
# New version for AGU - only plot VS
vel_i=1
for view_i in range(len(view)):
    rt.plot_rays(home,run_name,residname,veltype[vel_i],view[view_i],axlims[view_i],stations,events,by_path,mymap,faultfile)
    print 'Plotted velocity type '+str(veltype[vel_i])+', view '+str(view[view_i])
    plt.close()
    plt.close()
    
#Plot with the cutoff value:
#Plot in a loop
vel_i=1
for view_i in range(len(view)):
    rt.plot_rays_cutoffval(home,run_name,residname,veltype[vel_i],view[view_i],axlims[view_i],stations,events,mymap,faultfile,cutoff_val)
    print 'Plotted cuttoff value '+str(cutoff_val)+', velocity type '+str(veltype[vel_i])+', view '+str(view[view_i])
    plt.close()
    plt.close()
        
#For the map view, plot the 3d raypaths once, then rotate and save:
axlims_3d=[[-116.9,-116.35],[33.3,33.75],[-20,1]]
vtype=2

##plot:
#figure3d=rt.plot_3d_raypaths(home,run_name,vtype,stations,events,axlims_3d,mymap,faultfile)


###################################
#########Grid path term############
#
#rpath='/media/vsahakian/katmai/anza/models/residuals/abdb_5sta_0-6.5_topography/abdb_5sta_0-6.5_topography_robj_raydat_interp_indices_FangVs.pckl'
#
##Info from propgrid (in the future write a function to do this automatically):
#nz=36 
#nlat=218 
#nlon=278
##bins:
#b=(nlon,nlat,nz)
#
#stattype='mean'
#veltype=1
#
#
#gobjpath=home+run_name+'/'+'abdb_5sta_0-6.5_topography_pterm_grid.pckl'
#gfile=open(gobjpath,'r')
#gobj=pickle.load(gfile)
#gfile.close()


