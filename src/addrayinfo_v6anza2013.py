import run_res
from os import path
import dread
import cPickle as pickle
import raytracing as rt
import matplotlib.pyplot as plt

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


#############################################################################
### FILES ###
#############

# V2anza2013 raytracing output files:
v2vs_rayfile = HOME+'/anza/fm3d/run_raytracing_v2anza2013_res4_Vs/rays.dat'
v2vp_rayfile = HOME+'/anza/fm3d/run_raytracing_v2anza2013_res4_Vp/rays.dat'

#v2anza residuals object:
v2path = ''
v3path = '/Users/vsahakian/anza/models/residuals/mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj.pckl'



#############################################################################


# Parse the rayfile:
print 'parsing vs for v2anza...'
v2_vs_path_list,v2_vs_receiver_id,v2_vs_source_id = rt.parse_rayfile(v2vs_rayfile)
print 'parsing vp for v2anza...'
v2_vp_path_list,v2_vp_receiver_id,v2_vp_source_id = rt.parse_rayfile(v2vp_rayfile)

# Open the residuals objects:
v2file = open(v2path,'r')
v2obj = pickle.load(v2file)
v2file.close()

v3file = open(v3path,'r')
v3obj = pickle.load(v3file)
v3file.close()