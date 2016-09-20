#Write out the locations of events and stations used in database
#to a file for plotting in GMT
#VJS 9/2016

import cPickle as pickle
from numpy import c_,savetxt

#What are the file locations to read/write?
dbpath='/Users/vsahakian/anza/data/abdb_5sta.pckl'

#GMT events:
evpath='/Users/vsahakian/anza/data/db_5sta_event_locations.ll'

#GMT stations:
stpath='/Users/vsahakian/anza/data/db_5sta_sta_locations.ll'

####
#Read in the database:
dbfile=open(dbpath,'r')
db=pickle.load(dbfile)
dbfile.close()

#Get output:
event_output=c_[db.elon, db.elat]
station_output=c_[db.stlon, db.stlat]


###
#Save output:

#format...
output_fmt='%9.5f\t%7.5f'

#save events
savetxt(evpath,event_output,output_fmt)
#and stations
savetxt(stpath,station_output,output_fmt)