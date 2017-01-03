##  Prune 2013 data....  ##

import cPickle as pickle
from numpy import where,unique, zeros
import cdefs as cdf
import dread as dr

## Read in file:
robjpath = '/media/vsahakian/katmai/anza/models/residuals/v2anza2013_Mc8.5_pgrid_5sta/v2anza2013_Mc8.5_pgrid_5sta_robj.pckl'

# Original database path name:
dbpath_in=''
# Out path:
dbpath_out=''

#  Find where residuals are greater than or less than a certain value:
resvalue = 4


###################3

#Open it into robj:
rfile=open(robjpath,'r')
robj=pickle.load(rfile)
rfile.close()

# Find where they're greater/less than value:
resless4=where(abs(robj.total_residual)<=4)[0]
resgreater4=where(abs(robj.total_residual)>4)[0]

# Find the length of this, and what percentage of the original database would be removed:
reslen=len(resgreater4)

# Get all values of databse for these indices:
# The basics...
evnum=robj.evnum[resgreater4]
elat=robj.elat[resgreater4]
elon=robj.elon[resgreater4]
edepth=robj.edepth[resgreater4]
sta=robj.sta[resgreater4]
stnum=robj.stnum[resgreater4]
ml=robj.ml[resgreater4]
mw=robj.mw[resgreater4]
pga=robj.pga[resgreater4]
pgv=robj.pgv[resgreater4]
pga_pg=robj.pga_pg[resgreater4]
r=robj.r[resgreater4]
vs30=robj.vs30[resgreater4]
ffdf=robj.ffdf[resgreater4]
md_ffdf=robj.md_ffdf[resgreater4]
stlat=robj.stlat[resgreater4]
stlon=robj.stlon[resgreater4]
stelv=robj.stelv[resgreater4]
        
# The residuals...
total_residual=robj.total_residual[resgreater4]
E_residual=robj.E_residual[resgreater4]
E_mean=robj.E_mean[resgreater4]
E_std=robj.E_std[resgreater4]
W_residual=robj.W_residual[resgreater4]
W_mean=robj.W_mean[resgreater4]
W_std=robj.W_std[resgreater4]
site_terms=robj.site_terms[resgreater4]
site_mean=robj.site_mean[resgreater4]
site_std=robj.site_std[resgreater4]
path_terms=robj.path_terms[resgreater4]
path_mean=robj.path_mean[resgreater4]
path_std=robj.path_std[resgreater4]
    

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


# Sample the original database
dr.recording_sample(dbpath_in,resless4,dbpath_out)

# Make the rejected residuals into a residual object of its own:
init_type='notbasic'
new_robj = cdf.residuals(None,None,None,init_style=init_type,evnum=evnum,elat=elat,elon=elon,edepth=edepth,sta=sta,stnum=stnum,ml=ml,mw=mw,
                    pga=pga,pgv=pgv,pga_pg=pga_pg,r=r,vs30=vs30,ffdf=ffdf,md_ffdf=md_ffdf,stlat=stlat,
                    stlon=stlon,stelv=stelv,source_i=source_i,receiver_i=receiver_i,total_residual=total_residual,E_residual=E_residual,
                    E_mean=E_mean,E_std=E_std,W_residual=W_residual,W_mean=W_mean,W_std=W_std,site_terms=site_terms,site_mean=site_mean,site_std=site_std,
                    path_terms=path_terms,path_mean=path_mean,path_std=path_std)

# Save that residual object:
