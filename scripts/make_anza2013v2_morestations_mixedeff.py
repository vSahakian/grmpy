###  Make database with more stations than events  ###
###  Sample by requiring each event to be recorded on some very large number of stations ###
###  VJS 1/2017

## Input Info
what_home=0

if what_home==0:
    #Desktop:
    HOME='/media/vsahakian/katmai'
elif what_home==1:
    #Mac:
    HOME='/Users/vsahakian'

num_stations = 31

dbpath_in = HOME+'/anza/data/databases/v2anza2013/v2anza2013_pgrid_5sta_res4.pckl'
dbpath_out = HOME+'/anza/data/databases/v2anza2013/v2anza2013_pgrid_'+str(num_stations)+'sta_res4.pckl'

#########################

import dread as dr
import cPickle as pickle
from numpy import unique

# Sample
dr.db_station_sample(dbpath_in,num_stations,dbpath_out)


# Read out back in and make some figures:
dbfile=open(dbpath_out,'r')
dbout=pickle.load(dbfile)
dbfile.close()

print 'There are now %i events, and %i stations' % (len(unique(dbout.evnum)),len(unique(dbout.sta)))

