##Get info for Alan Yong

import cPickle as pickle
from numpy import unique,savetxt,array,c_


#########   A   ########
###Get the info for the old database and stations
##sta lon lat elev vs30 site term


rpath_old='/media/vsahakian/katmai/anza/models/residuals/abdb_5sta_0-6.5_topography/abdb_5sta_0-6.5_topography_robj.pckl'
vs30_v1_ofile='/media/vsahakian/katmai/anza/data/vs30/vs30_database1.txt'

#Open object:
rfile=open(rpath_old,'r')
robjv1=pickle.load(rfile)
rfile.close()

#Get info:
sta_array,sta_ind=unique(robjv1.sta,return_index=True)
#stations are an array of arrays, for some reason...make it a list:
#zero otu list:
sta=[]
for station in range(len(sta_array)):
    sta.append(sta_array[station][0])

#Get other info...
lon=robjv1.stlon[sta_ind]
lat=robjv1.stlat[sta_ind]
elv=robjv1.stelv[sta_ind]
vs30=robjv1.vs30[sta_ind]
#Site terms is for some reason a list, not an array, so convert...
sitet=array(robjv1.site_terms)[sta_ind]

#concatenate and save:
outinfo=c_[sta,lon,lat,elv,vs30,sitet]
format='%s\t%s\t%s\t%s\t%s\t%s'
head='sta    lon    lat    elv    vs30    siteterm'
savetxt(vs30_v1_ofile,outinfo,fmt=format,header=head)


#########   B   ########
###Get the info for the new database and stations
##sta lon lat elev vs30 

rpath_new='/media/vsahakian/katmai/anza/data/databases/db2013_test/db2013test_5sta.pckl'
vs30_v2_ofile='/media/vsahakian/katmai/anza/data/vs30/vs30_database2.txt'

#Open
rfile=open(rpath_new,'r')
robj=pickle.load(rfile)
rfile.close()

#Get info:
sta,sta_ind=unique(robj.sta,return_index=True)

#Get other info...
lon=robj.stlon[sta_ind]
lat=robj.stlat[sta_ind]
elv=robj.stelv[sta_ind]
vs30=robj.vs30[sta_ind]

#concatenate and save:
outinfo=c_[sta,lon,lat,elv,vs30]
format='%s\t%s\t%s\t%s\t%s'
head='sta        lon        lat        elv         vs30    '
savetxt(vs30_v2_ofile,outinfo,fmt=format,header=head)