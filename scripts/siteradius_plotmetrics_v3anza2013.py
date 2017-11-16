###########################################################################
## Plot metrics for set distances away from teh site
## Compare them with path terms

import res_analysis as ra
import matplotlib.pyplot as plt
import pandas as pd
import cPickle as pickle
import numpy as np


what_home=0

if what_home==0:
    #Desktop:
    HOME='/katmai'
elif what_home==1:
    #Mac:
    HOME='/Users/vsahakian'
    
    
########################
#######  Paths  ########
########################

home=HOME+'/anza/models/residuals/'

metricdfpath = home + 'mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30_siteradiusmetrics.pckl'



########################
#######  Params  #######
########################

site_radius = np.array([0.5,1,2,5,10,20])
metric_list = ['path','normpath','gradpath']
colorscheme = 'viridis'

binplot_dims = [48, 16]
binplot_font = 14

###############################
######  Read in models  #######
###############################

dffile = open(metricdfpath,'r')
metricdf = pickle.load(dffile)
dffile.close()


###############################
########     Plots     ########
###############################

# First do for site term, by distance:
residualterm = 'site'
bin_by = 'Rrup'
axlims = [[0, 7], [-3, 5.5]]
color_by = 'Rrup'
clims = [0, 140, 1]
cbartick = 40
plotrowscols = [2, 4]


for radiusi in range(len(site_radius)):
    for metrici in range(len(metric_list)):
        metric = metric_list[metrici] + np.str(site_radius[radiusi])
        binedges = np.array([  0,  20,  40,  60,  80, 100, 120, 180])
        
        i_binned_figure = ra.plot_binned_metric_dataframe(metricdf,residualterm,metric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,binplot_dims,plotrowscols,binplot_font)