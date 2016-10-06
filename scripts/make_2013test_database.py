#Make a database with Janine's test 2013 data:

import dread as dr
import cdefs as cdf
import numpy as np
import inversion as inv
import pickle
import matplotlib.pyplot as plt
import gmpe as gm


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

############
############

#Janine's flatfile:
flatfile=HOME+'/anza/data/databases/PGA_2013.dat'
#Vs30 file:
vs30file=HOME+'/anza/data/vs30/California.xyz'

#Output Vs30 file:
vs30name='vs30_2013test.txt'
vs30_outfile=HOME+'/anza/data/vs30/'+vs30name

#Read in the data:
evnum,evlat,evlon,evdep,sta,stlat,stlon,stelv,grcircle,ml,mw,pga_mgal=dr.read_jsbfile(flatfile)

#Compute Rrup for the data:
rrup=dr.compute_rrup(evlon,evlat,evdep,stlon,stlat,stelv)

######################
#Get Vs30:
#Find unique stations so it can be saved to a file later:
unique_stations=np.unique(sta)
unique_station_index=np.where(sta==unique_stations)

#Get the lat and lon of these points only:
stlon_unique=stlon[unique_station_index]
stlat_unique=stlat[unique_station_index]

vs30=dr.interp_vs30(stlat_unique,stlon_unique,vs30file)

#Save them to a file:
save_vs30=np.c_[unique_stations,stlon_unique,stlat_unique,vs30]
vs30fmt='%.4f\t%.4f\t%.4f'
np.savetxt(vs30_outfile,save_vs30,fmt=vs30fmt)

#Now read in the vs30 for each recording:

#######################


