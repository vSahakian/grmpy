##Play around iwth the raytracing...
#VJS 8/2016

import run_res
from os import path
import dread
import cPickle as pickle
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

home=HOME+'/anza/models/residuals/'
run_name='abdb_0-6.5_addindex'
dbpath=HOME+'/anza/data/abdb.pckl'
modelpath=HOME+'/anza/models/pckl/regr_0.0_1.0_2.0_6.5_VR_98.9.pckl'

#Set up the sources.in and receivers.in files for Vs.  So far, Vp has been run:
#in /media/vsahakian/katmai/anza/data/vm/fulltest_Vp

#Run these with raytracing write_sourcein and write_receiverin
#Use the directories defined in lines 15 - 33
#Import raytracing above
#Set velocity type:
veltype=2

#Write source.in file:
rt.write_sourcein(home,run_name,veltype)

#Write recievers.in file:
rt.write_receiverin(home,run_name)


######RUN##########

#For vp:
rayfile='/media/vsahakian/katmai/anza/data/vm/fulltest_Vp/rays.dat'
veltype=1

#Read in:
rt.store_rayinfo(home,run_name,rayfile,veltype)

####

#For vs:
rayfile='/media/vsahakian/katmai/anza/data/vm/fulltest_Vs/rays.dat'
veltype=2

#Read in:
rt.store_rayinfo(home,run_name,rayfile,veltype)

###
#Plot rays:
def plot_rays(home,run_name,veltype,view,axlims,stations,events)