##  Prune 2013 data....  ##

import cPickle as pickle
from numpy import where,unique, zeros

## Read in file:
robjpath = '/media/vsahakian/katmai/anza/models/residuals/v2anza2013_Mc8.5_pgrid_5sta/v2anza2013_Mc8.5_pgrid_5sta_robj.pckl'

#  Find where residuals are greater than or less than a certain value:
resvalue = 4


###################3

#Open it into robj:
rfile=open(robjpath,'r')
robj=pickle.load(rfile)
rfile.close()

# Find where they're greater/less than value:
reswhere = where(abs(robj.total_residual)>4)[0]

# Find the length of this, and what percentage of the original database would be removed:
reslen = len(reswhere)

# Get all values of databse for these indices:
# The basics...
evnum=robj.evnum[reswhere]
elat=robj.elat[reswhere]
elon=robj.elon[reswhere]
edepth=robj.edepth[reswhere]
sta=robj.sta[reswhere]
stnum=robj.stnum[reswhere]
ml=robj.ml[reswhere]
mw=robj.mw[reswhere]
pga=robj.pga[reswhere]
pgv=robj.pgv[reswhere]
pga_pg=robj.pga_pg[reswhere]
r=robj.r[reswhere]
vs30=robj.vs30[reswhere]
ffdf=robj.ffdf[reswhere]
md_ffdf=robj.md_ffdf[reswhere]
stlat=robj.stlat[reswhere]
stlon=robj.stlon[reswhere]
stelv=robj.stelv[reswhere]
        
# The residuals...
total_residual=robj.total_residual[reswhere]
E_residual=robj.E_residual[reswhere]
E_mean=robj.E_mean[reswhere]
E_std=robj.E_std[reswhere]
W_residual=robj.W_residual[reswhere]
W_mean=robj.W_mean[reswhere]
W_std=robj.W_std[reswhere]
site_terms=robj.site_terms[reswhere]
site_mean=robj.site_mean[reswhere]
site_std=robj.site_std[reswhere]
path_terms=robj.path_terms[reswhere]
path_mean=robj.path_mean[reswhere]
path_std=robj.path_std[reswhere]
    

# Redo source and receiver......
# Get the unique station and event indices:
unique_events=unique(evnum)

# Zero out source ind array:
source_ind=zeros((len(evnum)))
# For each event in the record, devent, give it the source index to be used:
for event_ind in range(len(unique_events)):
    eventi=unique_events[event_ind]
    
    # Find where in the recordings list the event number is the same as this one:
    recording_event_ind=where(evnum==eventi)
    
    # Set the source ind to be one plus this event, so it indexes with the raytracing program:
    source_ind[recording_event_ind]=event_ind+1
    
# Now set these to integers...
source_i=source_ind.astype('int64')

##
# Next stations:
unique_stations=unique(stnum)

# Zero out array:
receiver_ind=zeros((len(stnum)))
# Loop through the unique stations:
for station_ind in range(len(unique_stations)):
    stationi=unique_stations[station_ind]
    
    # Find where in the recordings list the station is the same as this one:
    recording_station_ind=where(stnum==stationi)[0]

    # Set the receiver ind to be one plus this station, so it indexes with the raytracin gprogram:
    receiver_ind[recording_station_ind]=station_ind+1    
    
# Set these to integers:
receiver_i=receiver_ind.astype('int64')


# Print the  number of events and stations involved in this portion of hte data:
nsta_old = len(unique(robj.sta))
nsta_new = len(unique(sta))

nev_old = len(unique(robj.evnum))
nev_new = len(unique(evnum))

nrecords_old = len(robj.total_residual)
nrecords_new = len(total_residual)
percent_db_kept = ((nrecords_old - nrecords_new)/nrecords_old)*100

print 'There were %i stations; after pruning there are %i stations' % (nsta_old,nsta_new)
print 'There were %i events; after pruning there are %i events' % (nev_old,nev_new)
print 'Removing residuals greater than abs(%f) leaves %f percent of the data' % (resvalue,percent_db_kept)

# Make it into a residual object of its own:


# Save that residual object:
   