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

#home=HOME+'/anza/models/residuals/'
#run_name='abdb_0-6.5_addindex'
#dbpath=HOME+'/anza/data/abdb.pckl'
#modelpath=HOME+'/anza/models/pckl/regr_0.0_1.0_2.0_6.5_VR_98.9.pckl'

home=HOME+'/anza/models/residuals/'
run_name='abdb_5sta_0-6.5_VR'
dbpath=HOME+'/anza/data/abdb_5sta.pckl'
modelpath=HOME+'/anza/models/pckl/regr_0.0_6.5_VR_98.9.pckl'

########################
####Set up Raytracing###
########################

##Set up the sources.in and receivers.in files for Vs.  So far, Vp has been run:
##in /media/vsahakian/katmai/anza/data/vm/fulltest_Vp
#
##Run these with raytracing write_sourcein and write_receiverin
##Use the directories defined in lines 15 - 33
##Import raytracing above
##Set velocity type:
#veltype=2
#
##Write source.in file:
#rt.write_sourcein(home,run_name,veltype)
#
##Write recievers.in file:
#rt.write_receiverin(home,run_name)


#####################
###Store ray files###
#####################

#For vp:
rayfile=HOME+'/anza/data/vm/fulltest_Vp/rays.dat'
veltype=1

#Read in:
rt.store_rayinfo(home,run_name,rayfile,veltype)

#####

#For vs:
rayfile=HOME+'/anza/data/vm/fulltest_Vs/rays.dat'
veltype=2

#Read in:
rt.store_rayinfo(home,run_name,rayfile,veltype)

###
#Plot rays:
#Vp:
veltype=[1,2]
#map:
view=[0,1,2]
axlims=[[-116.9,-116.35],[33.3,33.75]]
stations=1
events=1
by_path=1
mymap='jet'


###
#Plot in a loop
for vel_i in range(len(veltype)):
    for view_i in range(len(view)):
        rt.plot_rays(home,run_name,veltype[vel_i],view[view_i],axlims,stations,events,by_path,mymap)
        plt.close()
        plt.close()
