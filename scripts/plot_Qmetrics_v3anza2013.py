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
run_name='mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30/'


## Path to residuals object:
rpath=HOME+'/anza/models/residuals/mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30/mixedcoeff_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_robj_FILTERED_raydat_interp_indices_HaukQp_Qs.pckl'

## Figure directory:
fig_dir = HOME+'/anza/models/residuals/mixedregr_v3anza2013_pga_5coeff_a4_-1.20_Mc_8.5_res4_noVs30/figs/HaukQ'


###############################
##### Read in objects #########
###############################

## Read in the residuals object
rfile = open(rpath,'r')
robj = pickle.load(rfile)
rfile.close()


#Get some correlation coefficeints:

path_qp_pears_coeff,ptail=pearsonr(robj.path_terms,robj.ind_s_qp_pathint)
normpath_qp_pears_coeff,nptail=pearsonr(robj.path_terms,robj.ind_s_qp_normpathint)
grad_pears_qp_coeff,gptail=pearsonr(robj.path_terms,robj.ind_s_qp_gradpathint)

path_qs_pears_coeff,ptail=pearsonr(robj.path_terms,robj.ind_s_qs_pathint)
normpath_qs_pears_coeff,nptail=pearsonr(robj.path_terms,robj.ind_s_qs_normpathint)
grad_pears_qs_coeff,gptail=pearsonr(robj.path_terms,robj.ind_s_qs_gradpathint)



########################
#######  Plots  ########
########################

#Some more plots...
index=['ind_s_qp_pathint','ind_s_qp_normpathint','ind_s_qp_gradpathint','ind_s_qs_pathint','ind_s_qs_normpathint','ind_s_qs_gradpathint']
color_by=['r','mw']
axlims=[[[-4,4.5],[420,220000]],[[-4,4.5],[0.01,0.16]],[[-4,4.5],[10,910]],[[-4,4.5],[500,333000]],[[-4,4.5],[0.01,0.17]],[[-4,4.5],[10,820]]]
cmap='viridis'
cvals=[[0,150],[0.5,3.2]]

for indexi in range(len(index)):
    for color_by_i in range(len(color_by)):
        ra.plot_pathterms_colored(home,run_name,robj,index[indexi],axlims[indexi],color_by[color_by_i],cvals[color_by_i],cmap,force_fig_dir=fig_dir)

#Now for other terms...
term='site_terms'
index=['ind_s_qp_pathint','ind_s_qp_normpathint','ind_s_qp_gradpathint','ind_s_qs_pathint','ind_s_qs_normpathint','ind_s_qs_gradpathint']
color_by=['r','mw']
axlims=[[[-4,4.5],[420,220000]],[[-4,4.5],[0.01,0.16]],[[-4,4.5],[10,910]],[[-4,4.5],[500,333000]],[[-4,4.5],[0.01,0.17]],[[-4,4.5],[10,820]]]
cmap='viridis'
cvals=[[0,150],[0.5,3.2]]

for indexi in range(len(index)):
    for color_by_i in range(len(color_by)):
        ra.plot_terms_colored(home,run_name,robj,term,index[indexi],axlims[indexi],color_by[color_by_i],cvals[color_by_i],cmap,force_fig_dir=fig_dir)
        plt.close('all')

#Event residual:
term='E_residual'
index=['ind_s_qp_pathint','ind_s_qp_normpathint','ind_s_qp_gradpathint','ind_s_qs_pathint','ind_s_qs_normpathint','ind_s_qs_gradpathint']
color_by=['r','mw']
axlims=[[[-4,4.5],[420,220000]],[[-4,4.5],[0.01,0.16]],[[-4,4.5],[10,910]],[[-4,4.5],[500,333000]],[[-4,4.5],[0.01,0.17]],[[-4,4.5],[10,820]]]
cmap='viridis'
cvals=[[0,150],[0.5,3.2]]

for indexi in range(len(index)):
    for color_by_i in range(len(color_by)):
        ra.plot_terms_colored(home,run_name,robj,term,index[indexi],axlims[indexi],color_by[color_by_i],cvals[color_by_i],cmap,force_fig_dir=fig_dir)
        plt.close('all')



##########################################################################
######################       Path      ###################################
##########################################################################

## Make binned distance plot for Qs and gradient metric

# Bin and color gradient/path term by distance:
plotterm = 'path'
plotmetric = 'ind_s_qs_gradpathint'
binedges = array([  0,  20,  40,  60,  80, 100, 120, 180])
bin_by = 'Rrup'
axlims = [[10, 815], [-3, 5.5]]
color_by = 'Rrup'
colorscheme = 'viridis'
clims = [0, 140, 1]
cbartick = 40
plotdims = [48, 16]
plotrowscols = [2, 4]
fontsz = 14

figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by

path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz)

path_distance_distance.savefig(fig_dir + '/pdfs/' + figname + '.pdf')
path_distance_distance.savefig(fig_dir + '/' + figname + '.png')


# Bin and color gradient/path term by magnitude:
plotterm = 'path'
plotmetric = 'ind_s_qs_gradpathint'
binedges = array([  0, 1.0,  1.5,  2.0, 2.5, 3.0, 4.0])
bin_by = 'M'
axlims = [[10, 815], [-3, 5.5]]
color_by = 'M'
colorscheme = 'viridis'
clims = [0, 3.05, 0.1]
cbartick = 1.0
plotdims = [48, 16]
plotrowscols = [2, 3]
fontsz = 14

figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by

path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz)

path_distance_distance.savefig(fig_dir + '/pdfs/' + figname + '.pdf')
path_distance_distance.savefig(fig_dir + '/' + figname + '.png')


# Bin gradient/path term by distance, color by M:
plotterm = 'path'
plotmetric = 'ind_s_qs_gradpathint'
binedges = array([  0,  20,  40,  60,  80, 100, 120, 180])
bin_by = 'Rrup'
axlims = [[10, 815], [-3, 5.5]]
color_by = 'M'
colorscheme = 'viridis'
clims = [0, 3.05, 0.1]
cbartick = 1
plotdims = [48, 16]
plotrowscols = [2, 4]
fontsz = 14

figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by

path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz)

path_distance_distance.savefig(fig_dir + '/pdfs/' + figname + '.pdf')
path_distance_distance.savefig(fig_dir + '/' + figname + '.png')


#############################################################################

## Make plots for Qp and gradient metric

# Bin and color gradient/path term by distance:
plotterm = 'path'
plotmetric = 'ind_s_qp_gradpathint'
binedges = array([  0,  20,  40,  60,  80, 100, 120, 180])
bin_by = 'Rrup'
axlims = [[10, 920], [-3, 5.5]]
color_by = 'Rrup'
colorscheme = 'viridis'
clims = [0, 140, 1]
cbartick = 40
plotdims = [48, 16]
plotrowscols = [2, 4]
fontsz = 14

figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by

path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz)

path_distance_distance.savefig(fig_dir + '/pdfs/' + figname + '.pdf')
path_distance_distance.savefig(fig_dir + '/' + figname + '.png')


# Bin and color gradient/path term by magnitude:
plotterm = 'path'
plotmetric = 'ind_s_qp_gradpathint'
binedges = array([  0, 1.0,  1.5,  2.0, 2.5, 3.0, 4.0])
bin_by = 'M'
axlims = [[10, 920], [-3, 5.5]]
color_by = 'M'
colorscheme = 'viridis'
clims = [0, 3.05, 0.1]
cbartick = 1.0
plotdims = [48, 16]
plotrowscols = [2, 3]
fontsz = 14

figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by

path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz)

path_distance_distance.savefig(fig_dir + '/pdfs/' + figname + '.pdf')
path_distance_distance.savefig(fig_dir + '/' + figname + '.png')

# Bin gradient/path term by distance, color by M:
plotterm = 'path'
plotmetric = 'ind_s_qp_gradpathint'
binedges = array([  0,  20,  40,  60,  80, 100, 120, 180])
bin_by = 'Rrup'
axlims = [[10, 920], [-3, 5.5]]
color_by = 'M'
colorscheme = 'viridis'
clims = [0, 3.05, 0.1]
cbartick = 1
plotdims = [48, 16]
plotrowscols = [2, 4]
fontsz = 14

figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by

path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz)

path_distance_distance.savefig(fig_dir + '/pdfs/' + figname + '.pdf')
path_distance_distance.savefig(fig_dir + '/' + figname + '.png')



##########################################################################

## Make binned distance plot for Qs and norm path metric

# Bin and color gradient/path term by distance:
plotterm = 'path'
plotmetric = 'ind_s_qs_normpathint'
binedges = array([  0,  20,  40,  60,  80, 100, 120, 180])
bin_by = 'Rrup'
axlims = [[0.04, 0.17], [-3, 5.5]]
color_by = 'Rrup'
colorscheme = 'viridis'
clims = [0, 140, 1]
cbartick = 40
plotdims = [48, 16]
plotrowscols = [2, 4]
fontsz = 26
plt_xticks = array([0.04,0.08,0.12,0.16])

figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by

path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz,xticks=plt_xticks)

path_distance_distance.savefig(fig_dir + '/pdfs/' + figname + '.pdf',transparent=True)
path_distance_distance.savefig(fig_dir + '/' + figname + '.png',transparent=True)


# Bin and color gradient/path term by magnitude:
plotterm = 'path'
plotmetric = 'ind_s_qs_normpathint'
binedges = array([  0, 1.0,  1.5,  2.0, 2.5, 3.0, 4.0])
bin_by = 'M'
axlims = [[0.04, 0.17], [-3, 5.5]]
color_by = 'M'
colorscheme = 'viridis'
clims = [0, 3.05, 0.1]
cbartick = 1.0
plotdims = [48, 16]
plotrowscols = [2, 3]
fontsz = 26
plt_xticks = array([0.04,0.08,0.12,0.16])

figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by

path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz,xticks=plt_xticks)

path_distance_distance.savefig(fig_dir + '/pdfs/' + figname + '.pdf',transparent=True)
path_distance_distance.savefig(fig_dir + '/' + figname + '.png',transparent=True)


# Bin gradient/path term by distance, color by M:
plotterm = 'path'
plotmetric = 'ind_s_qs_normpathint'
binedges = array([  0,  20,  40,  60,  80, 100, 120, 180])
bin_by = 'Rrup'
axlims = [[0.04, 0.17], [-3, 5.5]]
color_by = 'M'
colorscheme = 'viridis'
clims = [0, 3.05, 0.1]
cbartick = 1
plotdims = [48, 16]
plotrowscols = [2, 4]
fontsz = 26
plt_xticks = array([0.04,0.08,0.12,0.16])

figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by

path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz,xticks=plt_xticks)

path_distance_distance.savefig(fig_dir + '/pdfs/' + figname + '.pdf')
path_distance_distance.savefig(fig_dir + '/' + figname + '.png')


#############################################################################

## Make binned distance plot for Qp and norm path metric

# Bin and color gradient/path term by distance:
plotterm = 'path'
plotmetric = 'ind_s_qp_normpathint'
binedges = array([  0,  20,  40,  60,  80, 100, 120, 180])
bin_by = 'Rrup'
axlims = [[0.0, 0.17], [-3, 5.5]]
color_by = 'Rrup'
colorscheme = 'viridis'
clims = [0, 140, 1]
cbartick = 40
plotdims = [48, 16]
plotrowscols = [2, 4]
fontsz = 26
plt_xticks = array([0.04,0.08,0.12,0.16])

figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by

path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz,xticks=plt_xticks)

path_distance_distance.savefig(fig_dir + '/pdfs/' + figname + '.pdf',transparent=True)
path_distance_distance.savefig(fig_dir + '/' + figname + '.png',transparent=True)


# Bin and color gradient/path term by magnitude:
plotterm = 'path'
plotmetric = 'ind_s_qp_normpathint'
binedges = array([  0, 1.0,  1.5,  2.0, 2.5, 3.0, 4.0])
bin_by = 'M'
axlims = [[0.0, 0.17], [-3, 5.5]]
color_by = 'M'
colorscheme = 'viridis'
clims = [0, 3.05, 0.1]
cbartick = 1.0
plotdims = [48, 16]
plotrowscols = [2, 3]
fontsz = 26
plt_xticks = array([0.04,0.08,0.12,0.16])

figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by

path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz,xticks=plt_xticks)

path_distance_distance.savefig(fig_dir + '/pdfs/' + figname + '.pdf',transparent=True)
path_distance_distance.savefig(fig_dir + '/' + figname + '.png',transparent=True)


# Bin gradient/path term by distance, color by M:
plotterm = 'path'
plotmetric = 'ind_s_qp_normpathint'
binedges = array([  0,  20,  40,  60,  80, 100, 120, 180])
bin_by = 'Rrup'
axlims = [[0.00, 0.17], [-3, 5.5]]
color_by = 'M'
colorscheme = 'viridis'
clims = [0, 3.05, 0.1]
cbartick = 1
plotdims = [48, 16]
plotrowscols = [2, 4]
fontsz = 26
plt_xticks = array([0.04,0.08,0.12,0.16])

figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by

path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz,xticks=plt_xticks)

path_distance_distance.savefig(fig_dir + '/pdfs/' + figname + '.pdf')
path_distance_distance.savefig(fig_dir + '/' + figname + '.png')




#############################################################################
#
#
## Bin gradient/path term by distance, color by azimuth:
#plotterm = 'path'
#plotmetric = 'ind_s_qs_gradpathint'
#binedges = array([  0,  20,  40,  60,  80, 100, 120, 180])
#bin_by = 'Rrup'
#axlims = [[10, 815], [-3, 5.5]]
#color_by = 'Site Azimuth'
#colorscheme = 'viridis'
#clims = [0, 360, 1]
#cbartick = 90
#plotdims = [48, 16]
#plotrowscols = [2, 4]
#fontsz = 14
#
#figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by
#
#path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz)
#
#path_distance_distance.savefig(fig_dir + 'pdfs/' + figname + '.pdf')
#path_distance_distance.savefig(fig_dir + figname + '.png')
#
#
## Bin gradient/path term by azimuth, color by distance:
#plotterm = 'path'
#plotmetric = 'ind_s_vs_gradpathint'
#binedges = array([  0,  45,  90,  135,  180,  225, 270, 315, 360])
#bin_by = 'Site Azimuth'
#axlims = [[0, 7], [-3, 5.5]]
#color_by = 'Rrup'
#colorscheme = 'viridis'
#clims = [0, 140, 1]
#cbartick = 40
#plotdims = [48, 16]
#plotrowscols = [2, 4]
#fontsz = 14
#
#figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by
#
#path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz)
#
#path_distance_distance.savefig(fig_dir + 'pdfs/' + figname + '.pdf')
#path_distance_distance.savefig(fig_dir + figname + '.png')
#
#
##################
#
#
## Bin and color gradient/path term normalized by distance, by distance:
#plotterm = 'path'
#plotmetric = 'ind_s_vs_gradpathint_normdist'
#binedges = array([  0,  20,  40,  60,  80, 100, 120, 180])
#bin_by = 'Rrup'
#axlims = [[0, 7], [-3, 5.5]]
#color_by = 'Rrup'
#colorscheme = 'viridis'
#clims = [0, 140, 1]
#cbartick = 40
#plotdims = [48, 16]
#plotrowscols = [2, 4]
#fontsz = 14
#
#figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by
#
#path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz)
#
#path_distance_distance.savefig(fig_dir + 'pdfs/' + figname + '.pdf')
#path_distance_distance.savefig(fig_dir + figname + '.png')




##########################################################################
######################       Site      ###################################
##########################################################################

# Plot site stuff:

# Bin gradient/site term by distance, color by distance:
plotterm = 'site'
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

figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by

path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz)

path_distance_distance.savefig(fig_dir + 'pdfs/' + figname + '.pdf')
path_distance_distance.savefig(fig_dir + figname + '.png')



# Bin distance-normalized gradient/site term by distance, color by distance:
plotterm = 'site'
plotmetric = 'ind_s_vs_gradpathint_normdist'
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

figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by

path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz)

path_distance_distance.savefig(fig_dir + 'pdfs/' + figname + '.pdf')
path_distance_distance.savefig(fig_dir + figname + '.png')




# Bin gradient/site term by distance, color by distance with 10km spacing:
plotterm = 'site'
plotmetric = 'ind_s_vs_gradpathint'
binedges = array([  0,  10,  20])
bin_by = 'Rrup'
axlims = [[0, 7], [-3, 5.5]]
color_by = 'Rrup'
colorscheme = 'viridis'
clims = [0, 19, 1]
cbartick = 4
plotdims = [48, 16]
plotrowscols = [1, 2]
fontsz = 14

figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by + '_10kmspacing'

path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz)

path_distance_distance.savefig(fig_dir + 'pdfs/' + figname + '.pdf')
path_distance_distance.savefig(fig_dir + figname + '.png')





# Bin gradient/site term by azimuth, color by azimuth:
plotterm = 'site'
plotmetric = 'ind_s_vs_gradpathint'
binedges = array([  0,  45,  90,  135,  180,  225, 270, 315, 360])
bin_by = 'Site Azimuth'
axlims = [[0, 7], [-3, 5.5]]
color_by = 'Site Azimuth'
colorscheme = 'viridis'
clims = [0, 360, 1]
cbartick = 90
plotdims = [48, 16]
plotrowscols = [2, 4]
fontsz = 14

figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by

path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz)

path_distance_distance.savefig(fig_dir + 'pdfs/' + figname + '.pdf')
path_distance_distance.savefig(fig_dir + figname + '.png')


# Bin gradient/site term by distance, color by azimuth:
plotterm = 'site'
plotmetric = 'ind_s_vs_gradpathint'
binedges = array([  0,  20,  40,  60,  80, 100, 120, 180])
bin_by = 'Rrup'
axlims = [[0, 7], [-3, 5.5]]
color_by = 'Site Azimuth'
colorscheme = 'viridis'
clims = [0, 360, 1]
cbartick = 90
plotdims = [48, 16]
plotrowscols = [2, 4]
fontsz = 14

figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by

path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz)

path_distance_distance.savefig(fig_dir + 'pdfs/' + figname + '.pdf')
path_distance_distance.savefig(fig_dir + figname + '.png')



# Bin gradient/site term by distance, color by magnitude:
plotterm = 'site'
plotmetric = 'ind_s_vs_gradpathint'
binedges = array([  0,  20,  40,  60,  80, 100, 120, 180])
bin_by = 'Rrup'
axlims = [[0, 7], [-3, 5.5]]
color_by = 'M'
colorscheme = 'viridis'
clims = [0, 3.05, 0.1]
cbartick = 1.0
plotdims = [48, 16]
plotrowscols = [2, 4]
fontsz = 14

figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by

path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz)

path_distance_distance.savefig(fig_dir + 'pdfs/' + figname + '.pdf')
path_distance_distance.savefig(fig_dir + figname + '.png')



########################################################################
### Site median/mean plots

# Bin mean gradient/site term by distance, color by distance
plotterm = 'site'
plotmetric = 'ind_s_vs_gradpathint_sitemean'
binedges = array([  0,  20,  40,  60,  80, 100, 120, 180])
bin_by = 'Rrup'
axlims = [[0, 7], [-3, 5.5]]
color_by = 'Rrup'
color_by_stat = 'mean'
colorscheme = 'viridis'
clims = [0, 140, 1]
cbartick = 40
plotdims = [48, 16]
plotrowscols = [2, 4]
fontsz = 14

figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by

path_distance_distance = ra.plot_binned_metric_statistic(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,color_by_stat,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz)

path_distance_distance.savefig(fig_dir + 'pdfs/' + figname + '.pdf')
path_distance_distance.savefig(fig_dir + figname + '.png')


# Bin median gradient/site term by distance, color by distance
plotterm = 'site'
plotmetric = 'ind_s_vs_gradpathint_sitemedian'
binedges = array([  0,  20,  40,  60,  80, 100, 120, 180])
bin_by = 'Rrup'
axlims = [[0, 7], [-3, 5.5]]
color_by = 'Rrup'
color_by_stat = 'median'
colorscheme = 'viridis'
clims = [0, 140, 1]
cbartick = 40
plotdims = [48, 16]
plotrowscols = [2, 4]
fontsz = 14

figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by

path_distance_distance = ra.plot_binned_metric_statistic(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,color_by_stat,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz)

path_distance_distance.savefig(fig_dir + 'pdfs/' + figname + '.pdf')
path_distance_distance.savefig(fig_dir + figname + '.png')






##########################################################################
######################   Path and Site ###################################
##########################################################################

# Plot path and site stuff:

# Bin gradient/path and site term by distance, color by distance:
plotterm = 'path_site'
plotmetric = 'ind_s_vs_gradpathint'
binedges = array([  0,  20,  40,  60,  80, 100, 120, 180])
bin_by = 'Rrup'
axlims = [[0, 7], [-3.5, 5.5]]
color_by = 'Rrup'
colorscheme = 'viridis'
clims = [0, 140, 1]
cbartick = 40
plotdims = [48, 16]
plotrowscols = [2, 4]
fontsz = 14

figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by

path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz)

path_distance_distance.savefig(fig_dir + 'pdfs/' + figname + '.pdf')
path_distance_distance.savefig(fig_dir + figname + '.png')



# Bin gradient/path and site term by distance, color by azimuth:
plotterm = 'path_site'
plotmetric = 'ind_s_vs_gradpathint'
binedges = array([  0,  20,  40,  60,  80, 100, 120, 180])
bin_by = 'Rrup'
axlims = [[0, 7], [-3.5, 5.5]]
color_by = 'Site Azimuth'
colorscheme = 'viridis'
clims = [0, 360, 1]
cbartick = 90
plotdims = [48, 16]
plotrowscols = [2, 4]
fontsz = 14

figname = plotmetric + '_' + plotterm + '_' + bin_by + '_' + color_by

path_distance_distance = ra.plot_binned_metric(robj,plotterm,plotmetric,binedges,bin_by,axlims,color_by,colorscheme,clims,cbartick,plotdims,plotrowscols,fontsz)

path_distance_distance.savefig(fig_dir + 'pdfs/' + figname + '.pdf')
path_distance_distance.savefig(fig_dir + figname + '.png')
