## Make metric plots
# VJS 10/2017


import dread
import cdefs as cdf
import cPickle as pickle
from numpy import meshgrid,zeros,where,array
import run_res_analysis as runra
import res_analysis as ra
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr

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

########################
#######  Paths  ########
########################

home=HOME+'/anza/models/residuals/'
run_name='mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30'


## Path to residuals object:
rpath=HOME+'/anza/models/residuals/mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_raydat_interp_indices_FangVs.pckl'


###############################
##### Read in objects #########
###############################

## Read in the residuals object
rfile = open(rpath,'r')
robj = pickle.load(rfile)
rfile.close()


#Get some correlation coefficeints:
path_pears_coeff,ptail=pearsonr(robj.path_terms,robj.ind_s_vs_pathint)
normpath_pears_coeff,nptail=pearsonr(robj.path_terms,robj.ind_s_vs_normpathint)
grad_pears_coeff,gptail=pearsonr(robj.path_terms,robj.ind_s_vs_gradpathint)



########################
#######  Plots  ########
########################

#Some more plots...
index=['ind_s_vs_pathint','ind_s_vs_normpathint','ind_s_vs_gradpathint']
color_by=['r','mw']
axlims=[[[-4,4.5],[10,8450]],[[-4,4.5],[0.60,0.85]],[[-4,4.5],[0,6.5]]]
cmap='viridis'
cvals=[[0,150],[0.5,3.2]]

for indexi in range(len(index)):
    for color_by_i in range(len(color_by)):
        ra.plot_pathterms_colored(home,run_name,robj,index[indexi],axlims[indexi],color_by[color_by_i],cvals[color_by_i],cmap)

#Now for other terms...
term='site_terms'
index=['ind_s_vs_pathint','ind_s_vs_normpathint','ind_s_vs_gradpathint']
color_by=['r','mw']
axlims=[[[-4,4.5],[10,8450]],[[-4,4.5],[0.60,0.85]],[[-4,4.5],[0,6.5]]]
cmap='viridis'
cvals=[[0,150],[0.5,3.2]]

for indexi in range(len(index)):
    for color_by_i in range(len(color_by)):
        ra.plot_terms_colored(home,run_name,robj,term,index[indexi],axlims[indexi],color_by[color_by_i],cvals[color_by_i],cmap)
        

#Event residual:
term='E_residual'
index=['ind_s_vs_pathint','ind_s_vs_normpathint','ind_s_vs_gradpathint']
color_by=['r','mw']
axlims=[[[-4,4.5],[10,8450]],[[-4,4.5],[0.60,0.85]],[[-4,4.5],[0,6.5]]]
cmap='viridis'
cvals=[[0,150],[0.5,3.2]]

for indexi in range(len(index)):
    for color_by_i in range(len(color_by)):
        ra.plot_terms_colored(home,run_name,robj,term,index[indexi],axlims[indexi],color_by[color_by_i],cvals[color_by_i],cmap)




########
## Make binned distance plot for gradient metric

# Bin and color by distance:
plotterm = 'path'
plotmetric = 'ind_s_vs_gradpathint'
binedges = array([  0,  20,  40,  60,  80, 100, 120, 180])
bin_by = 'Rrup'
axlims = [[0, 7], [-3, 5.5]]
color_by = 'Rrup'
colorscheme = 'viridis'
clims = [0, 140, 1]
cbartick = 40
plotdims = [48, 16]
plotrowscols = [2, 4]
fontsz = 14

path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz)

