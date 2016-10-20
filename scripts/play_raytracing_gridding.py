##Play around iwth the raytracing...
#VJS 8/2016

import run_res
from os import path
import dread
import cPickle as pickle
import raytracing as rt
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import res_analysis as ra

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
run_name='abdb_5sta_0-6.5_topography'
dbpath=HOME+'/anza/data/abdb_5sta.pckl'
modelpath=HOME+'/anza/models/pckl/regr_0.0_6.5_VR_98.9.pckl'
faultfile=HOME+'/anza/data/faults/Holocene_LatestPleistocene_117.5w_115.5w_33n_34n.pckl'
materialpath=HOME+'/anza/data/pckl/FangVs.pckl'
gobjpath=home+run_name+'/'+'abdb_5sta_0-6.5_topography_pterm_grid.pckl'

##################################
########Grid path term############

#Info from propgrid (in the future write a function to do this automatically):
#bins:
nz=36 
nlat=218 
nlon=278
bindims=(nlon,nlat,nz)
#stat type, etc:
stattype='mean'
veltype=1

#Plotting parameters:
slicecoord=33.5
coordtype='lat'
cmap='jet'
aspectr=0.008
climits=[-1.2,1.2]
#slice_axlims=[[33.35,33.71],[-19.67,2.62]]
slice_axlims=[[-116.9,-116.4],[-19.67,2.62]]

#############################
#######Make grid object:#####
#############################

#gobj=ra.grid_path_term(home,run_name,bindims,veltype,stattype)

#############################
#############################

#############################
###...Or Open grid object:###
#############################

##Open the grid object:
gfile=open(gobjpath,'r')
gobj=pickle.load(gfile)
gfile.close()

#Open the materials object:
mfile=open(materialpath,'r')
mobj=pickle.load(mfile)
mfile.close()

###Make plots...
#INitiate plot
fig,axarr=plt.subplots(2,1)

#First plot the velocity model:


gridax=gobj.plot_slice(axarr[1],slicecoord,coordtype,aspectr,cmap,climits,slice_axlims)