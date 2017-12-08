###########################################################################
## Plot metrics for set distances away from teh site
## Compare them with path terms

import res_analysis as ra
import matplotlib.pyplot as plt
import pandas as pd
import cPickle as pickle
import numpy as np


what_home=1

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

metricdfpath = home + 'mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30/mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30_siteradiusmetrics.pckl'
metriccsvpath = home + 'mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30/mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30_siteradiusmetrics.csv'


fig_dir = home + 'mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30/figs/'
#fig_dir = '/Users/vsahakian/Desktop/debugfigs/'


########################
#######  Params  #######
########################

site_radius = np.array([0.5,1,2,5,10,20])
#site_radius = np.array([0.5])

# metric_list = ['path','normpath','gradpath']
metric_list = ['gradpath']
colorscheme = 'viridis'

binplot_font = 14
markersize = 30

###############################
######  Read in models  #######
###############################

#dffile = open(metricdfpath,'r')
#metricdf = pickle.load(dffile)
#dffile.close()

metricdf = pd.read_csv(metriccsvpath)


###############################
########     Plots     ########
###############################
#
## First do for site term, for one set distance bin:
#residualterm = ['site','path']
#bin_by = ['Rrup','Site Azimuth']
#binedges = [np.array([0,180]),np.array([0,360])]
##axlims = [[[10, 1240], [-3, 2.5]], [[0.5, 1], [-3, 2.5]], [[0, 5], [-3, 2.5]] ]
#axlims = [[[0, 5], [-3, 2.5]] ]
#
#color_by = ['Rrup','Site Azimuth']
#clims = [[0, 140, 1], [0, 360, 1]]
#cbartick = [40,90]
#plotrowscols = [1, 1]
#binplot_dims = [15,10]
#
#print 'Starting plotting...'
#for residuali in range(len(residualterm)):
#    print 'residual ' + residualterm[residuali]
#    
#    for radiusi in range(len(site_radius)):
#        print 'radius ' + np.str(site_radius[radiusi])
#        
#        for metrici in range(len(metric_list)):
#            print 'metric' + metric_list[metrici]
#            
#            for binbyi in range(len(bin_by)):
#                print 'binning by ' + bin_by[binbyi]
#
#                    
#                iresidual = residualterm[residuali]
#                ibinedges = binedges[binbyi]
#                ibin_by = bin_by[binbyi]
#                icolor_by = color_by[binbyi]
#                iclims = clims[binbyi]
#                icbartick = cbartick[binbyi]
#                imetric = metric_list[metrici] + np.str(site_radius[radiusi])
#                
#                i_binned_figure = ra.plot_binned_metric_dataframe(metricdf,iresidual,imetric,ibinedges,ibin_by,axlims[metrici],icolor_by,colorscheme,iclims,icbartick,binplot_dims,plotrowscols,binplot_font)
#                
#                i_figname = 'siteradius_' + iresidual + '_' + imetric + '_binby' + ibin_by + '_colorby' + icolor_by + '_' + np.str(len(ibinedges)-1) + 'bins'
#                i_figpath = fig_dir + i_figname + '.png'
#                i_pdfpath = fig_dir + 'pdfs/' + i_figname + '.pdf'
#        
#                i_binned_figure.savefig(i_figpath)
#                i_binned_figure.savefig(i_pdfpath)
#                
#                plt.close(i_binned_figure)


#######################################################################################

#residualterm = ['site']
#bin_by = ['Rrup','Site Azimuth']
#binedges = [np.array([0,5,10,20,40,60,180]),np.array([-0.00001,45,90,135,180,225,270,315,360])]
#plotrowscols = [[3,2],[4,2]]
##axlims = [[[10, 1240], [-3, 2.5]], [[0.5, 1], [-3, 2.5]], [[0, 5], [-3, 2.5]] ]
#axlims = [[[0, 4.5], [-3, 2.5]] ]
##axlims = [[[0, 1], [-3, 2.5]], [[0, 1], [-3, 2.5]], [[0, 1], [-3, 2.5]], [[0, 1], [-3, 2.5]], [[0, 1], [-3, 2.5]], [[0, 1], [-3, 2.5]]]
#
#color_by = ['Rrup','Site Azimuth']
#clims = [[0, 80, 1], [0, 360, 1]]
#cbartick = [40,90]
#binplot_dims = [[12,9.75],[12,13]]
#statistic = 'mean'
#
#for residuali in range(len(residualterm)):
#    for metrici in range(len(metric_list)):
#        for radiusi in range(len(site_radius)):
#            for binbyi in range(len(bin_by)):
#                i_residual = residualterm[residuali]
#                i_metric = metric_list[metrici] + np.str(site_radius[radiusi])
#                i_axlims = axlims[metrici]
#                i_bin_by = bin_by[binbyi]
#                i_binedges = binedges[binbyi]
#                i_plotrc = plotrowscols[binbyi]
#                i_plotdims = binplot_dims[binbyi]
#                i_color_by = color_by[binbyi]
#                i_clims = clims[binbyi]
#                i_cbartick = cbartick[binbyi]
#                
#                print i_residual
#                print i_metric
#                print i_bin_by
#                
#                
#                i_statfig = ra.plot_binned_metric_statistic_df(metricdf,i_residual,i_metric,statistic,i_binedges,i_bin_by,i_axlims,i_color_by,statistic,colorscheme,i_clims,i_cbartick,i_plotdims,i_plotrc,binplot_font,markersize)
#                
#                if i_bin_by == 'Site Azimuth':
#                    binby_text = 'SiteAzimuth'
#                else:
#                    binby_text = i_bin_by
#                    
#                if i_color_by == 'Site Azimuth':
#                    colorby_text = 'Site Azimuth'
#                else:
#                    colorby_text = i_color_by
#                    
#                i_figname = 'siteradius_' + i_residual + '_' + i_metric + '_binby' + binby_text + '_colorby' + colorby_text + '_' + statistic + '_' + np.str(len(i_binedges)-1) + 'bins'
#                i_figpath = fig_dir + i_figname + '.png'
#                i_pdfpath = fig_dir + 'pdfs/' + i_figname + '.pdf'
#        
#                i_statfig.savefig(i_figpath)
#                i_statfig.savefig(i_pdfpath)
#                
#                plt.close(i_statfig)                
                
                
######################################################################################
# First do for site term, for one set distance bin:
residualterm = ['site','path']
bin_by = ['Rrup','Site Azimuth']
binedges = [np.array([0,180]),np.array([0,360])]
#axlims = [[[10, 1240], [-3, 2.5]], [[0.5, 1], [-3, 2.5]], [[0, 5], [-3, 2.5]] ]
axlims = [[[0, 5], [-3, 2.5]] ]

color_by = ['Rrup','Site Azimuth']
clims = [[0, 140, 1], [0, 360, 1]]
cbartick = [40,90]
plotrowscols = [1, 1]
binplot_dims = [15,10]

print 'Starting plotting...'
for residuali in range(len(residualterm)):
    print 'residual ' + residualterm[residuali]
    
    for radiusi in range(len(site_radius)):
        print 'radius ' + np.str(site_radius[radiusi])
        
        for metrici in range(len(metric_list)):
            print 'metric' + metric_list[metrici]
            
            for binbyi in range(len(bin_by)):
                print 'binning by ' + bin_by[binbyi]

                    
                iresidual = residualterm[residuali]
                ibinedges = binedges[binbyi]
                ibin_by = bin_by[binbyi]
                icolor_by = color_by[binbyi]
                iclims = clims[binbyi]
                icbartick = cbartick[binbyi]
                imetric = metric_list[metrici] + np.str(site_radius[radiusi])
                
                i_binned_figure = ra.plot_binned_metric_dataframe(metricdf,iresidual,imetric,ibinedges,ibin_by,axlims[metrici],icolor_by,colorscheme,iclims,icbartick,binplot_dims,plotrowscols,binplot_font)
                
                i_figname = 'siteradius_' + iresidual + '_' + imetric + '_binby' + ibin_by + '_colorby' + icolor_by + '_' + np.str(len(ibinedges)-1) + 'bins'
                i_figpath = fig_dir + i_figname + '.png'
                i_pdfpath = fig_dir + 'pdfs/' + i_figname + '.pdf'
        
                i_binned_figure.savefig(i_figpath)
                i_binned_figure.savefig(i_pdfpath)
                
                plt.close(i_binned_figure)
