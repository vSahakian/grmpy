import run_res
from os import path
import dread
import cPickle as pickle
import raytracing as rt
import matplotlib.pyplot as plt
import numpy as np

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
v2path = '/Users/vsahakian/anza/models/residuals/mixedregr_v2anza2013_pga_5coeff_Mc_8.5_res4_noVs30/mixedcoeff_v2anza2013_pga_noVs30_5coeff_pga__ncoeff5_Mc_8.5_VR_99.9_robj.pckl'
v3path = '/Users/vsahakian/anza/models/residuals/mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj.pckl'



#############################################################################

############################
####### Open things ########
############################

# Parse the rayfile:
print 'parsing vs for v2anza...'
v2_vs_path_list,v2_vs_receiver_id,v2_vs_source_id = rt.parse_rayfile(v2vs_rayfile)
print 'parsing vp for v2anza...'
v2_vp_path_list,v2_vp_receiver_id,v2_vp_source_id = rt.parse_rayfile(v2vp_rayfile)

# Open the residuals objects:
print 'opening v2anza residuals object'
v2file = open(v2path,'r')
v2obj = pickle.load(v2file)
v2file.close()

print 'opening v3anza residuals object'
v3file = open(v3path,'r')
v3obj = pickle.load(v3file)
v3file.close()

############################

## REad in the v3 database. For every recording in it, find that recording inside of 
#   v2, and then get that source and receiver id. Look for that combo in the path 
#   list, and read in that ray data.

# First initiate the path lists:
v3_vp_pathlist = []
v3_vs_pathlist = []

# Also initiate a list to add to if which v3 database recordings have a raypath:
v3_withraypath_ind = []

for i_v3recording in range(len(v3obj.pga_pg)):
    i_ev_v3 = v3obj.evnum[i_v3recording]
    i_st_v3 = v3obj.sta[i_v3recording]
    
    # Find where this exists in v2:
    v3record_in_v2_ind = np.where((v2obj.evnum == i_ev_v3) & (v2obj.sta == i_st_v3))[0]
    
    # If this recording has a matching recording in database v2, keep it and the rest of its info:
    if len(v3record_in_v2_ind)>0:
        # Append to list:
        v3_withraypath_ind.append(i_v3recording)
    
        # Get source and receiver numbers for this:
        i_v2_sourcei = v2obj.source_i[v3record_in_v2_ind]
        i_v2_receiveri = v2obj.receiver_i[v3record_in_v2_ind]
    
        # Find where these are inside the path list, and grab those:
        v3record_in_pathlist_ind = np.where((v2_vs_receiver_id == i_v2_receiveri) & (v2_vs_source_id == i_v2_sourcei))[0][0]
    
        # Now get the path list:
        v3_vp_pathlist.append(v2_vp_path_list[v3record_in_pathlist_ind])
        v3_vs_pathlist.append(v2_vs_path_list[v3record_in_pathlist_ind])

    
### NEXT::
# 1. convertr lon, lat, depth to appropriate measures with right lon_type
# 2. Make new residuals object with the recordings that have raypaths
# 3. Add the raypaths to it
# 4. Save it....come up with a useful name...
