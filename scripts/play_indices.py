##
#Wrapper script for material model stuff...

import dread
import cdefs as cdf
import cPickle as pickle
from numpy import meshgrid,zeros,where
import run_res_analysis as runra
import res_analysis as ra
from matplotlib.pyplot import scatter

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









########
########












###############
##############

##Read in residuals object
#rfile=open(rpath,'r')
#robj=pickle.load(rfile)
#rfile.close()
#
##Read in data:
#x,y,z,nx,ny,nz,mod=dread.read_material_model(coordspath,materialmodelpath)
#
##Convert longitude:
#x=x-360
#
##Make object:
#mmodel=cdf.material_model(x,y,z,nx,ny,nz,mod)
#
#
#####Test for the first ray:#######
##Get actual values of x, y, and z:
#grid_x=mmodel.x
#grid_y=mmodel.y
#grid_z=-mmodel.z
#
##Get the values to interpolate (aka values of ray...for now just one ray):
#ray_x=robj.vs_lon[0]
#ray_y=robj.vs_lat[0]
#ray_z=robj.vs_depth[0]
#
#
###
##Make points for griddata:
#gX,gY,gZ=meshgrid(grid_x,grid_y,grid_z,indexing='ij')
#
##Make column vectors
#columnx=gX.ravel()
#columny=gY.ravel()
#columnz=gZ.ravel()
#data=mmodel.materials.transpose(2,1,0).ravel()
#
##Make interpolator
#interpolator = Rbf(columnx, columny, columnz, data,function='linear')
#
#######
###Interpolate
#######
#
##First initialize empty list to put interpolated values in:
#ray_dat=[]
#
#for ray_i in range(len(robj.vs_lon)):
#    ray_x_i=robj.vs_lon[ray_i]
#    ray_y_i=robj.vs_lat[ray_i]
#    ray_z_i=robj.vs_depth[ray_i]
#    
#    #Interpolate:
#    vs_i=interpolator(ray_x_i,ray_y_i,ray_z_i)
#    
#    #Append to the list:
#    ray_dat.append(vs_i)
#    
########
###Interpolate by reshaping:
###
#
#
#
#######
###Single ray:
#ray_x=robj.vs_lon[0]
#ray_y=robj.vs_lat[0]
#ray_z=robj.vs_depth[0]
#
#vs=interpolator(ray_x,ray_y,ray_z)
#
###Scatter ray:
#i=where(columnz==-3)[0]
##Scatter interpolater velocity model:
#scatter(columnx[i],columny[i],c=data[i],s=70)
#
##Scatter ray and velocity on ray:
#scatter(ray_x,ray_z,c=vs,s=60,lw=0.1)
#
