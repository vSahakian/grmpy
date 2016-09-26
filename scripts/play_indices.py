##
#Wrapper script for material model stuff...

import dread
import cdefs as cdf

#Paths to Fang 2016 model
coordspath='/media/vsahakian/katmai/anza/data/vm/Fang2016/coords.txt'
materialmodelpath='/media/vsahakian/katmai/anza/data/vm/Fang2016/Vs.dat'

#Read in data:
x,y,z,nx,ny,nz,mod=dread.read_material_model(coordspath,materialmodelpath)

#Convert longitude:
x=x-360

#Make object:
materialmodel=cdf.material_model(x,y,z,nx,ny,nz,mod)


####Test for the first ray:#######
#Get actual values of x, y, and z:
grid_x=