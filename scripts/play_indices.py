##
#Wrapper script for material model stuff...

import dread
import cdefs as cdf
import cPickle as pickle
from numpy import meshgrid,zeros
from scipy.interpolate import griddata,Rbf

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
    

#Paths to Fang 2016 model
coordspath=HOME+'/anza/data/vm/Fang2016/coords.txt'
materialmodelpath=HOME+'/anza/data/vm/Fang2016/Vs.dat'

#Path to residuals object:
rpath=HOME+'/anza/models/residuals/abdb_5sta_0-6.5_topography/abdb_5sta_0-6.5_topography_robj_raydat.pckl'

########
#######

#Read in residuals object
rfile=open(rpath,'r')
robj=pickle.load(rfile)
rfile.close()

#Read in data:
x,y,z,nx,ny,nz,mod=dread.read_material_model(coordspath,materialmodelpath)

#Convert longitude:
x=x-360

#Make object:
mmodel=cdf.material_model(x,y,z,nx,ny,nz,mod)


####Test for the first ray:#######
#Get actual values of x, y, and z:
grid_x=mmodel.x
grid_y=mmodel.y
grid_z=-mmodel.z

#Get the values to interpolate (aka values of ray...for now just one ray):
ray_x=robj.vs_lon[0]
ray_y=robj.vs_lat[0]
ray_z=robj.vs_depth[0]


##
#Make points for griddata:
gX,gY,gZ=meshgrid(grid_x,grid_y,grid_z,indexing='ij')

#Make column vectors
columnx=gX.ravel()
columny=gY.ravel()
columnz=gZ.ravel()
data=mmodel.materials.transpose(2,1,0).ravel()

#Make interpolator
interpolator = Rbf(columnx, columny, columnz, data,function='linear')

#Interpolate
vs=interpolator(ray_x,ray_y,ray_z)

##Scatter ray:
scatter(ray_x,ray_z,c=vs,s=60,lw=0.1)

#Scatter interpolater velocity model:
scatter(columnx[i],columny[i],c=data[i],s=70)