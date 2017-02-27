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
    
#Paths:
home=HOME+'/anza/models/residuals/'
run_name='mixedregr_anza2013_Mc_8.5'

#Paths to Fang 2016 model
coordspath=HOME+'/anza/data/vm/Fang2016/coords.txt'
materialmodelpath=HOME+'/anza/data/vm/Fang2016/Vs.dat'

#Path to residuals object:
rpath=HOME+'/anza/models/residuals/mixedregr_anza2013_Mc_8.5/mixedregr_anza2013_Mc_8.5_VR_99.9_robj_raydat.pckl'

#Path to store materialobject:
materialmodelname='FangVs'

#Name of the material model object location
matobjpath=HOME+'/anza/data/pckl/'+materialmodelname+'.pckl'

#Material flag - what kind is it?  0=Vp, 1=Vs, 2=Vp/Vs, 3=Qp, 4=Qs, 5=Qp/Qs
materialflag=1

#Residual Object paths:
#interpolated ray values
rpath_interp=rpath.split('.pckl')[0]+'_interp_'+materialmodelname+'.pckl'

#interpolated ray values and indices
rpath_indices=rpath.split('.pckl')[0]+'_interp_indices_'+materialmodelname+'.pckl'

##############################################################################

###
##Interpolate rays:
#Read in residuals object
print 'opening residuals object...'
rfile=open(rpath_indices,'r')
robj=pickle.load(rfile)
rfile.close()
print 'residuals object opened'


###############
##############
## Even more plots....
index='ind_s_vs_gradpathint'
term='path_terms'
color_by=['r','mw']
axlims=[[-3,6],[0,6.5]]
cmap='jet'
cvals=[[0,170],[1,3.5]]

# Separate dataset:
condition='gr2mw'
gr2ind=where(robj.mw>=2)[0]
rx2=robj.path_terms[gr2ind]
ry2=robj.ind_s_vs_gradpathint[gr2ind]

for color_by_i in range(len(color_by)):
    ra.plot_terms_colored_condition(home,run_name,condition,robj,rx2,ry2,term,index,axlims,color_by[color_by_i],gr2ind,cvals[color_by_i],cmap)

####
condition='gr2le2_3km'
gr2le23ind=where((robj.mw>=2) & (robj.mw<=2.3))[0]
rx223=robj.path_terms[gr2le23ind]
ry223=robj.ind_s_vs_gradpathint[gr2le23ind]

for color_by_i in range(len(color_by)):
    ra.plot_terms_colored_condition(home,run_name,condition,robj,rx223,ry223,term,index,axlims,color_by[color_by_i],gr2le23ind,cvals[color_by_i],cmap)
    
#####
####
condition='gr2_3le2_5km'
gr23le25ind=where((robj.mw>=2.3) & (robj.mw<=2.5))[0]
rx2325=robj.path_terms[gr23le25ind]
ry2325=robj.ind_s_vs_gradpathint[gr23le25ind]

for color_by_i in range(len(color_by)):
    ra.plot_terms_colored_condition(home,run_name,condition,robj,rx2325,ry2325,term,index,axlims,color_by[color_by_i],gr23le25ind,cvals[color_by_i],cmap)
    
####
###
condition='gr2_5le2_7km'
gr25le27ind=where((robj.mw>=2.5) & (robj.mw<=2.7))[0]
rx2527=robj.path_terms[gr25le27ind]
ry2527=robj.ind_s_vs_gradpathint[gr25le27ind]

for color_by_i in range(len(color_by)):
    ra.plot_terms_colored_condition(home,run_name,condition,robj,rx2527,ry2527,term,index,axlims,color_by[color_by_i],gr25le27ind,cvals[color_by_i],cmap)
    
###
###
condition='gr2_7le2_9km'
gr27le29ind=where((robj.mw>=2.7) & (robj.mw<=2.9))[0]
rx2729=robj.path_terms[gr27le29ind]
ry2729=robj.ind_s_vs_gradpathint[gr27le29ind]

for color_by_i in range(len(color_by)):
    ra.plot_terms_colored_condition(home,run_name,condition,robj,rx2729,ry2729,term,index,axlims,color_by[color_by_i],gr27le29ind,cvals[color_by_i],cmap)

###
###
condition='gr2_9le3_1km'
gr29le31ind=where((robj.mw>=2.9) & (robj.mw<=3.1))[0]
rx2931=robj.path_terms[gr29le31ind]
ry2931=robj.ind_s_vs_gradpathint[gr29le31ind]

for color_by_i in range(len(color_by)):
    ra.plot_terms_colored_condition(home,run_name,condition,robj,rx2931,ry2931,term,index,axlims,color_by[color_by_i],gr29le31ind,cvals[color_by_i],cmap)
   
    
###
###
condition='gr3_1le3_3km'
gr31le33ind=where((robj.mw>=3.1) & (robj.mw<=3.3))[0]
rx3133=robj.path_terms[gr31le33ind]
ry3133=robj.ind_s_vs_gradpathint[gr31le33ind]

for color_by_i in range(len(color_by)):
    ra.plot_terms_colored_condition(home,run_name,condition,robj,rx3133,ry3133,term,index,axlims,color_by[color_by_i],gr31le33ind,cvals[color_by_i],cmap)    
    
    
###
###
condition='gr3_3km'
gr33ind=where(robj.mw>=3.3)[0]
rx33=robj.path_terms[gr33ind]
ry33=robj.ind_s_vs_gradpathint[gr33ind]

for color_by_i in range(len(color_by)):
    ra.plot_terms_colored_condition(home,run_name,condition,robj,rx33,ry33,term,index,axlims,color_by[color_by_i],gr33ind,cvals[color_by_i],cmap)