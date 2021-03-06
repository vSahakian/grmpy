##
#Wrapper script for material model stuff...

import dread
import cdefs as cdf
import cPickle as pickle
from numpy import meshgrid,zeros,where
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
run_name='abdb_5sta_0-6.5_topography'

#Paths to Fang 2016 model
coordspath=HOME+'/anza/data/vm/Fang2016/coords.txt'
materialmodelpath=HOME+'/anza/data/vm/Fang2016/Vs.dat'

#Path to residuals object:
rpath=HOME+'/anza/models/residuals/abdb_5sta_0-6.5_topography/abdb_5sta_0-6.5_topography_robj_raydat.pckl'

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

########
########

###
##Interpolate rays:
#Read in residuals object
rfile=open(rpath,'r')
robj=pickle.load(rfile)
rfile.close()

#Make material model:
#Convert the longitudes from positive west to negative:
lonconvert=2
materialmodel=runra.make_material_object(coordspath,materialmodelpath,matobjpath,lonconvert)


#Interpolate Vp and Vs rays through this model:
interpolation_type='linear'
p_ray_data,s_ray_data=runra.interp_rays(robj,materialmodel,interpolation_type)

#Add to and save in a residuals object:
#values for p rays:
robj.add_material_values(p_ray_data,materialflag,0)
#Values for s rays:
robj.add_material_values(s_ray_data,materialflag,1)

rfile_interp=open(rpath_interp,'w')
pickle.dump(robj,rfile_interp)
rfile_interp.close()


###
##Compute indices and save to a new residuals object:
#First compute for the p-waves:
#No normalization:
ind_p_vs_path=ra.compute_pathintegral(p_ray_data,materialmodel,0)
#normalization:
ind_p_vs_normpath=ra.compute_pathintegral(p_ray_data,materialmodel,1)
#gradient:
ind_p_vs_gradpath=ra.compute_devpathintegral(p_ray_data,materialmodel,0)

#Now for s-wves:
ind_s_vs_path=ra.compute_pathintegral(s_ray_data,materialmodel,0)
#normalization:
ind_s_vs_normpath=ra.compute_pathintegral(s_ray_data,materialmodel,1)
#gradient:
ind_s_vs_gradpath=ra.compute_devpathintegral(s_ray_data,materialmodel,0)


##Save the indices:
#indext type =0, ray type is p=0,material type is vs=1) 
robj.add_indices(ind_p_vs_path,0,0,1)
#indextype =1,
robj.add_indices(ind_p_vs_normpath,1,0,1)
#indextype=2:
robj.add_indices(ind_p_vs_gradpath,2,0,1)

#indext type =0, ray type is s=1,material type is vs=1) 
robj.add_indices(ind_s_vs_path,0,1,1)
#indextype =1,
robj.add_indices(ind_s_vs_normpath,1,1,1)
#indextype=2:
robj.add_indices(ind_s_vs_gradpath,2,1,1)


##Save to object:
#
rfile_indices=open(rpath_indices,'w')
pickle.dump(robj,rfile_indices)
rfile_indices.close()

#Get some correlation coefficeints:
grad_pears_coeff,gptail=pearsonr(robj.path_terms,robj.ind_s_vs_gradpathint)
normpath_pears_coeff,nptail=pearsonr(robj.path_terms,robj.ind_s_vs_normpathint)
path_pears_coeff,ptail=pearsonr(robj.path_terms,robj.ind_s_vs_pathint)

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
axlims=[[[-4,4],[200,1400]],[[-4,4],[0.66,0.81]],[[-4,4],[0,5]]]
cmap='jet'
cvals=[[0,20],[1.2,3]]

for indexi in range(len(index)):
    for color_by_i in range(len(color_by)):
        ra.plot_pathterms_colored(home,run_name,robj,index[indexi],axlims[indexi],color_by[color_by_i],cvals[color_by_i],cmap)

#Now for other terms...
term='site_terms'
index=['ind_s_vs_pathint','ind_s_vs_normpathint','ind_s_vs_gradpathint']
color_by=['r','mw']
axlims=[[[-1.5,1.5],[200,1400]],[[-1.5,1.5],[0.66,0.81]],[[-1.5,1.5],[0,5]]]
cmap='jet'
cvals=[[0,20],[1.2,3]]

for indexi in range(len(index)):
    for color_by_i in range(len(color_by)):
        ra.plot_terms_colored(home,run_name,robj,term,index[indexi],axlims[indexi],color_by[color_by_i],cvals[color_by_i],cmap)
        

#Event residual:
term='E_residual'
index=['ind_s_vs_pathint','ind_s_vs_normpathint','ind_s_vs_gradpathint']
color_by=['r','mw']
axlims=[[[-1.5,1.5],[200,1400]],[[-1.5,1.5],[0.66,0.81]],[[-1.5,1.5],[0,5]]]
cmap='jet'
cvals=[[0,20],[1.2,3]]

for indexi in range(len(index)):
    for color_by_i in range(len(color_by)):
        ra.plot_terms_colored(home,run_name,robj,term,index[indexi],axlims[indexi],color_by[color_by_i],cvals[color_by_i],cmap)

###############
##############

