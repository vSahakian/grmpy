import run_res
from os import path
import dread
import cPickle as pickle
import raytracing as rt
import matplotlib.pyplot as plt
import numpy as np
import cdefs as cdf

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
v3_vp_raydepth = []
v3_vp_raylon = []
v3_vp_raylat = []

v3_vs_raydepth = []
v3_vs_raylon = []
v3_vs_raylat = []

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
        v3_vp_raydepth.append(v2_vp_path_list[v3record_in_pathlist_ind][:,0] - 6371)
        v3_vp_raylat.append(np.degrees(v2_vp_path_list[v3record_in_pathlist_ind][:,1]))
        v3_vp_raylon.append(np.degrees(v2_vp_path_list[v3record_in_pathlist_ind][:,2]) - 360)

        v3_vs_raydepth.append(v2_vs_path_list[v3record_in_pathlist_ind][:,0] - 6371)
        v3_vs_raylat.append(np.degrees(v2_vs_path_list[v3record_in_pathlist_ind][:,1]))
        v3_vs_raylon.append(np.degrees(v2_vs_path_list[v3record_in_pathlist_ind][:,2]) - 360)
        

# Make new residuals object with only recordings that have raypaths:
dbpath = 'none'
event_list_path = 'none'
station_list_path = 'none'

evnum_cut = v3obj.evnum[v3_withraypath_ind]
elat_cut = v3obj.elat[v3_withraypath_ind]
elon_cut = v3obj.elon[v3_withraypath_ind]
edepth_cut = v3obj.edepth[v3_withraypath_ind]
sta_cut = v3obj.sta[v3_withraypath_ind]
stnum_cut = v3obj.stnum[v3_withraypath_ind]
stlat_cut = v3obj.stlat[v3_withraypath_ind]
stlon_cut = v3obj.stlon[v3_withraypath_ind]
stelv_cut = v3obj.stelv[v3_withraypath_ind]
ml_cut = v3obj.ml[v3_withraypath_ind]
mw_cut = v3obj.mw[v3_withraypath_ind]
pga_cut = v3obj.pga[v3_withraypath_ind]
pgv_cut = None
pga_pg_cut = v3obj.pga_pg[v3_withraypath_ind]
pga_snr_cut = None
pgv_snr_cut = None
r_cut = v3obj.r[v3_withraypath_ind]
vs30_cut = v3obj.vs30[v3_withraypath_ind]
vs30_method_cut = v3obj.vs30_method[v3_withraypath_ind]
ffdf_cut = v3obj.ffdf[v3_withraypath_ind]
md_ffdf_cut = v3obj.md_ffdf[v3_withraypath_ind]

total_residual_cut = v3obj.total_residual[v3_withraypath_ind]
E_residual_cut = v3obj.E_residual[v3_withraypath_ind]
E_mean_cut = v3obj.E_mean
E_std_cut = v3obj.E_std
W_residual_cut = v3obj.W_residual[v3_withraypath_ind]
W_mean_cut = v3obj.W_mean
W_std_cut = v3obj.W_std
site_terms_cut = v3obj.site_terms[v3_withraypath_ind]
site_mean_cut = v3obj.site_mean
site_stderr_cut = v3obj.site_stderr[v3_withraypath_ind]
site_std_cut = v3obj.site_std
path_terms_cut = v3obj.path_terms[v3_withraypath_ind]
path_mean_cut = v3obj.path_mean
path_std_cut = v3obj.path_std

# Use old source and receiver numbres...in case they're useful...
source_i_cut = v3obj.source_i[v3_withraypath_ind]
receiver_i_cut = v3obj.receiver_i[v3_withraypath_ind]


v3_rayfiltered_obj = cdf.residuals(dbpath,event_list_path,station_list_path,init_style='notbasic',
                    evnum=evnum_cut,elat=elat_cut,elon=elon_cut,edepth=edepth_cut,sta=sta_cut,stnum=stnum_cut,ml=ml_cut,mw=mw_cut,
                    pga=pga_cut,pgv=pgv_cut,pga_pg=pga_pg_cut,pga_snr=pga_snr_cut,pgv_snr=pgv_snr_cut,r=r_cut,vs30=vs30_cut,ffdf=ffdf_cut,md_ffdf=md_ffdf_cut,stlat=stlat_cut,
                    stlon=stlon_cut,stelv=stelv_cut,source_i=source_i_cut,receiver_i=receiver_i_cut,vs30_method=vs30_method_cut,total_residual=total_residual_cut,E_residual=E_residual_cut,
                    E_mean=E_mean_cut,E_std=E_std_cut,W_residual=W_residual_cut,W_mean=W_mean_cut,W_std=W_std_cut,site_terms=site_terms_cut,site_mean=site_mean_cut,site_stderr=site_stderr_cut,site_std=site_std_cut,
                    path_terms=path_terms_cut,path_mean=path_mean_cut,path_std=path_std_cut)
    

    
### NEXT::

# 3. Add the raypaths to it
# 4. Save it....come up with a useful name...
