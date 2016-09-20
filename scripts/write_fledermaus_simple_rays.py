###
#Make fledermaus raypaths file for Debi
#Write out "simple" raypaths - just origin and endpoint, since they do not 
#turn much
#VJS 9/2016
#

import cPickle as pickle
from numpy import around

#Residual object pickle file with raypaths:
rpath='/media/vsahakian/katmai/anza/models/residuals/abdb_5sta_0-6.5_VR/abdb_5sta_0-6.5_VR_robj_raydat.pckl'

#Outpout fledermaus file, with the format:
#Format: lon_event lat_event dep_event lon_sta lat_sta dep_sta path_residual
fledermaus='/media/vsahakian/katmai/anza/models/residuals/abdb_5sta_0-6.5_VR/abdb_5sta_0-6.5_fledermaus_raypaths.txt'

#########

#Open the residuals object:
rfile=open(rpath,'r')
robj=pickle.load(rfile)
rfile.close()

##
#Open the fledermaus output file, print the header line:
headerline='lon_event  lat_event  dep_event  lon_station  lat_station  dep_station  path_residual\n'

f=open(fledermaus,'w')
f.write(headerline)

#
#Now for every recording in teh database, write the beginning and endpoints of the raypath (event/station), and the path term:

for recording_i in range(len(robj.evnum)):
    #Get event info - round and convert to string:
    ev_lon='%.4f' % robj.elon[recording_i]
    ev_lat='%.4f' % robj.elat[recording_i]
    ev_depth='%.4f' % robj.edepth[recording_i]
    
    #Get station info - set the station depth to 0:
    st_lon='%.4f' % robj.stlon[recording_i]
    st_lat='%.4f' % robj.stlat[recording_i]
    st_depth='%.4f' % 0.0000
    
    #Get path term:
    pathterm='%.5f' % robj.path_terms[recording_i]
    
    #Write out...
    dataline=ev_lon+'\t'+ev_lat+'\t'+ev_depth+'\t'+st_lon+'\t'+st_lat+'\t'+st_depth+'\t'+pathterm+'\n'
    f.write(dataline)
    
#Close the file:
f.close()
    
    