#Write out the locations of events and stations used in database
#to a file for plotting in GMT
#VJS 9/2016

import cPickle as pickle
from numpy import c_,savetxt,unique

##What are the file locations to read/write?
#dbpath='/Users/vsahakian/anza/data/abdb_5sta.pckl'
#
##GMT events:
#evpath='/Users/vsahakian/anza/data/db_5sta_event_locations.ll'
#
##GMT stations:
#stpath='/Users/vsahakian/anza/data/db_5sta_sta_locations.ll'


#What are the file locations to read/write?
#dbpath='/media/vsahakian/katmai/anza/data/databases/db2013_test/db2013test_5sta.pckl'
#dbpath='/media/vsahakian/katmai/anza/data/databases/anza2013/anza2013_5sta.pckl'
#dbpath='/Users/vsahakian/anza/data/databases/anza2013/anza2013_pgrid_5sta.pckl'
dbpath='/Users/vsahakian/anza/data/databases/v2anza2013/v2anza2013_pgrid_5sta_res4.pckl'



##GMT events:
#evpath='/Users/vsahakian/anza/data/db2013test_5sta.ll'
##GMT stations:
#stpath='/Users/vsahakian/anza/data/db2013test_5sta.ll'

##GMT events:
#evpath='/media/vsahakian/katmai/anza/data/databases/db2013_test/db2013test_5sta_events.ll'
##GMT stations:
#stpath='/media/vsahakian/katmai/anza/data/databases/db2013_test/db2013test_5sta_stations.ll'


##GMT events:
#evpath='/media/vsahakian/katmai/anza/data/databases/anza2013/anza2013_5sta_events.ll'
##GMT stations:
#stpath='/media/vsahakian/katmai/anza/data/databases/anza2013/anza2013_5sta_stations.ll'

##GMT events:
#evpath='/Users/vsahakian/anza/data/databases/anza2013/anza2013_pgrid_5sta_events.ll'
##GMT stations:
#stpath='/Users/vsahakian/anza/data/databases/anza2013/anza2013_pgrid_5sta_stations.ll'

#GMT events:
evpath='/Users/vsahakian/anza/data/databases/v2anza2013/v2anza2013_pgrid_5sta_res4_events.ll'
#GMT stations:
stpath='/Users/vsahakian/anza/data/databases/v2anza2013/v2anza2013_pgrid_5sta_res4_stations.ll'
stpath_names='/Users/vsahakian/anza/data/databases/v2anza2013/v2anza2013_pgrid_5sta_res4_stations_names.ll'

####
#Read in the database:
dbfile=open(dbpath,'r')
db=pickle.load(dbfile)
dbfile.close()

#Get output:
event_output=c_[db.elon, db.elat]
station_output=c_[db.stlon, db.stlat]

# Also, save one for Luke with the station name:
# unique indices...
usta,usta_ind = unique(db.sta,return_index=True)
station_names = db.sta[usta_ind]

###
#Save output:

#format...
output_fmt='%9.5f\t%7.5f'
output_name_fmt='%s\t%9.5f\t%7.5f'

#save events
savetxt(evpath,event_output,output_fmt)
#and stations
savetxt(stpath,station_output,output_fmt)

#and stations with names:
# Open file
f = open(stpath_names,'w')
for stationi in range(len(usta)):
    writeline='%5s \t %9.5f \t %7.5f \n' % (station_names[stationi],station_output[stationi][0],station_output[stationi][1])
    f.write(writeline)
f.close()
    
    
