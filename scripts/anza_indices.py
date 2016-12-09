##
#Wrapper script for material model stuff...

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

what_home=0

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

#########################################################################
#########################################################################

###
##Interpolate rays:
#Read in residuals object
print 'opening residuals object...'
rfile=open(rpath,'r')
robj=pickle.load(rfile)
rfile.close()
print 'residuals object opened'

#Make material model:
#Convert the longitudes from positive west to negative:
print 'Making material model:'
lonconvert=2
materialmodel=runra.make_material_object(coordspath,materialmodelpath,matobjpath,lonconvert)

print 'material model made. \n  Interpolating model...'

#Interpolate Vp and Vs rays through this model:
interpolation_type='linear'
#Run only for Vs
raytypes=array([1])
s_ray_data=runra.interp_rays(robj,materialmodel,interpolation_type,raytypes)

print 'Material model made.  \n adding to the residuals object...'

#Add to and save in a residuals object:
#values for p rays:
#robj.add_material_values(p_ray_data,materialflag,0)
#Values for s rays:
robj.add_material_values(s_ray_data,materialflag,1)

print 'Saving new res object...'
rfile_interp=open(rpath_interp,'w')
pickle.dump(robj,rfile_interp)
rfile_interp.close()


###
##Compute indices and save to a new residuals object:
##First compute for the p-waves:
##No normalization:
#ind_p_vs_path=ra.compute_pathintegral(p_ray_data,materialmodel,0)
##normalization:
#ind_p_vs_normpath=ra.compute_pathintegral(p_ray_data,materialmodel,1)
##gradient:
#ind_p_vs_gradpath=ra.compute_devpathintegral(p_ray_data,materialmodel,0)

print 'Computing indices for Vs...'
#Now for s-wves:
ind_s_vs_path=ra.compute_pathintegral(s_ray_data,materialmodel,0)
print 'vs path done'
#normalization:
ind_s_vs_normpath=ra.compute_pathintegral(s_ray_data,materialmodel,1)
print 'vs normalized path done'
#gradient:
ind_s_vs_gradpath=ra.compute_devpathintegral(s_ray_data,materialmodel,0)
print 'vs gradient along path done'

##Save the indices:
##indext type =0, ray type is p=0,material type is vs=1) 
#robj.add_indices(ind_p_vs_path,0,0,1)
##indextype =1,
#robj.add_indices(ind_p_vs_normpath,1,0,1)
##indextype=2:
#robj.add_indices(ind_p_vs_gradpath,2,0,1)

print 'Adding indices to object...'
#indext type =0, ray type is s=1,material type is vs=1) 
robj.add_indices(ind_s_vs_path,0,1,1)
print 'vs path added'
#indextype =1,
robj.add_indices(ind_s_vs_normpath,1,1,1)
print 'vs normalized path added'
#indextype=2:
robj.add_indices(ind_s_vs_gradpath,2,1,1)
print 'vs grad path added'


##Save to object:
#
print 'saving to object...'
rfile_indices=open(rpath_indices,'w')
pickle.dump(robj,rfile_indices)
rfile_indices.close()

print 'saved to object...'

#Get some correlation coefficeints:
grad_pears_coeff,gptail=pearsonr(robj.path_terms,robj.ind_s_vs_gradpathint)
normpath_pears_coeff,nptail=pearsonr(robj.path_terms,robj.ind_s_vs_normpathint)
path_pears_coeff,ptail=pearsonr(robj.path_terms,robj.ind_s_vs_pathint)
print 'got correlation coefficients'

########
########

##Make some plots...automate this later...##

#Plot of path term on x axis, gradient of path on y
plt.figure()
plt.scatter(robj.path_terms,robj.ind_s_vs_gradpathint,facecolors='none',edgecolors='blue')
plt.xlabel('Path Term (ln residual)')
plt.ylabel('Path integral of gradient of Vs')
plt.title('Pearsons coefficient = '+str(grad_pears_coeff))

#Plot of path term on x axis, normalized path integral on y
plt.figure()
plt.scatter(robj.path_terms,robj.ind_s_vs_normpathint,facecolors='none',edgecolors='blue')
plt.xlabel('Path Term (ln residual)')
plt.ylabel('Normalized integral of of Vs')
plt.title('Pearsons coefficient = '+str(normpath_pears_coeff))

#Plot of path term on x axis, path integral on y
plt.figure()
plt.scatter(robj.path_terms,robj.ind_s_vs_pathint,facecolors='none',edgecolors='blue')
plt.xlabel('Path Term (ln residual)')
plt.ylabel('Path integral of Vs')
plt.title('Pearsons coefficient = '+str(path_pears_coeff))



###
#Some more plots...
index=['ind_s_vs_pathint','ind_s_vs_normpathint','ind_s_vs_gradpathint']
color_by=['r','mw']
axlims=[[[-3,6],[100,8650]],[[-3,6],[0.66,0.85]],[[-3,6],[0,6.5]]]
cmap='jet'
cvals=[[0,140],[1.2,3]]

for indexi in range(len(index)):
    for color_by_i in range(len(color_by)):
        ra.plot_pathterms_colored(home,run_name,robj,index[indexi],axlims[indexi],color_by[color_by_i],cvals[color_by_i],cmap)

#Now for other terms...
term='site_terms'
index=['ind_s_vs_pathint','ind_s_vs_normpathint','ind_s_vs_gradpathint']
color_by=['r','mw']
axlims=[[[-2.5,2.5],[100,8650]],[[-2.5,2.5],[0.66,0.85]],[[-2.5,2.5],[0,6.5]]]
cmap='jet'
cvals=[[0,140],[1.2,3]]

for indexi in range(len(index)):
    for color_by_i in range(len(color_by)):
        ra.plot_terms_colored(home,run_name,robj,term,index[indexi],axlims[indexi],color_by[color_by_i],cvals[color_by_i],cmap)
        

#Event residual:
term='E_residual'
index=['ind_s_vs_pathint','ind_s_vs_normpathint','ind_s_vs_gradpathint']
color_by=['r','mw']
axlims=[[[-5,5],[100,8650]],[[-5,5],[0.66,0.85]],[[-5,5],[0,6.5]]]
cmap='jet'
cvals=[[0,140],[1.2,3]]

for indexi in range(len(index)):
    for color_by_i in range(len(color_by)):
        ra.plot_terms_colored(home,run_name,robj,term,index[indexi],axlims[indexi],color_by[color_by_i],cvals[color_by_i],cmap)

###############
##############
## Even more plots....
index='ind_s_vs_gradpathint'
term='path_terms'
color_by=['r','mw']
axlims=[[-3,6],[0,6.5]]
cmap='jet'
cvals=[[0,140],[1.2,3]]
condition='gr40km'

# Separate dataset:
gr40ind=where(robj.r>=40)[0]
rx40=robj.path_terms[gr40ind]
ry40=robj.ind_s_vs_gradpathint[gr40ind]

for color_by_i in range(len(color_by)):
    ra.plot_terms_colored_condition(home,run_name,condition,robj,rx40,ry40,term,index,axlims,color_by[color_by_i],gr40ind,cvals[color_by_i],cmap)

####
condition='gr40le60km'
gr40le60ind=where((robj.r>=40) & (robj.r<=60))[0]
rx4060=robj.path_terms[gr40le60ind]
ry4060=robj.ind_s_vs_gradpathint[gr40le60ind]

for color_by_i in range(len(color_by)):
    ra.plot_terms_colored_condition(home,run_name,condition,robj,rx4060,ry4060,term,index,axlims,color_by[color_by_i],gr40le60ind,cvals[color_by_i],cmap)
    
#####
####
condition='gr60le80km'
gr60le80ind=where((robj.r>=60) & (robj.r<=80))[0]
rx6080=robj.path_terms[gr60le80ind]
ry6080=robj.ind_s_vs_gradpathint[gr60le80ind]

for color_by_i in range(len(color_by)):
    ra.plot_terms_colored_condition(home,run_name,condition,robj,rx6080,ry6080,term,index,axlims,color_by[color_by_i],gr60le80ind,cvals[color_by_i],cmap)
    
####
###
condition='gr80le100km'
gr80le100ind=where((robj.r>=80) & (robj.r<=100))[0]
rx80100=robj.path_terms[gr80le100ind]
ry80100=robj.ind_s_vs_gradpathint[gr80le100ind]

for color_by_i in range(len(color_by)):
    ra.plot_terms_colored_condition(home,run_name,condition,robj,rx80100,ry80100,term,index,axlims,color_by[color_by_i],gr80le100ind,cvals[color_by_i],cmap)
    
###
###
condition='gr100le120km'
gr100le120ind=where((robj.r>=100) & (robj.r<=120))[0]
rx100120=robj.path_terms[gr100le120ind]
ry100120=robj.ind_s_vs_gradpathint[gr100le120ind]

for color_by_i in range(len(color_by)):
    ra.plot_terms_colored_condition(home,run_name,condition,robj,rx100120,ry100120,term,index,axlims,color_by[color_by_i],gr100le120ind,cvals[color_by_i],cmap)

###
###
condition='gr120km'
gr120ind=where(robj.r>=120)[0]
rx120=robj.path_terms[gr120ind]
ry120=robj.ind_s_vs_gradpathint[gr120ind]

for color_by_i in range(len(color_by)):
    ra.plot_terms_colored_condition(home,run_name,condition,robj,rx120,ry120,term,index,axlims,color_by[color_by_i],gr120ind,cvals[color_by_i],cmap)