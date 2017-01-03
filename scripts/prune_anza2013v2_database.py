##  Prune 2013 data....  ##

import cPickle as pickle
from numpy import where,unique, zeros,array,mean,std,setdiff1d
import cdefs as cdf
import dread as dr
import matplotlib.pyplot as plt

## Read in file:
robjpath = '/media/vsahakian/katmai/anza/models/residuals/v2anza2013_Mc8.5_pgrid_5sta/v2anza2013_Mc8.5_pgrid_5sta_robj.pckl'

# Output rejected residuals file:
robjpath_out = '/media/vsahakian/katmai/anza/models/residuals/v2anza2013_Mc8.5_pgrid_5sta/v2anza2013_Mc8.5_pgrid_5sta_robj_rejected.pckl'

# Original database path name:
dbpath_in='/media/vsahakian/katmai/anza/data/databases/v2anza2013/v2anza2013_pgrid_5sta.pckl'
# Out path:
dbpath_out='/media/vsahakian/katmai/anza/data/databases/v2anza2013/v2anza2013_pgrid_5sta_res4.pckl'

# Figure directory:
fig_dir='/media/vsahakian/katmai/anza/models/residuals/v2anza2013_Mc8.5_pgrid_5sta/figs/'

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
E_residual=array(robj.E_residual)[resgreater4]
E_mean=mean(E_residual)
E_std=std(E_residual)
W_residual=array(robj.W_residual)[resgreater4]
W_mean=mean(W_residual)
W_std=std(W_residual)
site_terms=array(robj.site_terms)[resgreater4]
site_mean=mean(site_terms)
site_std=std(site_terms)
path_terms=array(robj.path_terms)[resgreater4]
path_mean=mean(path_terms)
path_std=std(path_terms)
    

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
percent_db_kept = ((nrecords_old - nrecords_new)/float(nrecords_old))*100

print 'There were %i stations; after pruning there were rejections from %i unique stations' % (nsta_old,nsta_new)
print 'There were %i events; after pruning there were rejections from %i unique events' % (nev_old,nev_new)
print 'Removing residuals greater than abs(%f) leaves %f percent of the data' % (resvalue,percent_db_kept)


# Sample the original database
print 'Making sampled database'
dr.recording_sample(dbpath_in,resless4,dbpath_out)

#   Get and keep unique events and stations:
removed_events=setdiff1d(unique(robj.evnum),unique(robj.evnum[resless4]))
removed_stas=setdiff1d(unique(robj.sta),unique(robj.sta[resless4]))

# Make the rejected residuals into a residual object of its own:
init_type='notbasic'
new_robj = cdf.residuals(None,None,None,init_style=init_type,evnum=evnum,elat=elat,elon=elon,edepth=edepth,sta=sta,stnum=stnum,ml=ml,mw=mw,
                    pga=pga,pgv=pgv,pga_pg=pga_pg,r=r,vs30=vs30,ffdf=ffdf,md_ffdf=md_ffdf,stlat=stlat,
                    stlon=stlon,stelv=stelv,source_i=source_i,receiver_i=receiver_i,total_residual=total_residual,E_residual=E_residual,
                    E_mean=E_mean,E_std=E_std,W_residual=W_residual,W_mean=W_mean,W_std=W_std,site_terms=site_terms,site_mean=site_mean,site_std=site_std,
                    path_terms=path_terms,path_mean=path_mean,path_std=path_std)

# Save that residual object:
rout = open(robjpath_out,'w')
pickle.dump(new_robj,rout)
rout.close()


#######################
# Make some plots:

## EVENTS  ##
# Fig basename:
evname='rejected_event'

# First get the lone event(s) that were completely removed:
event_remove_ind=where(new_robj.evnum == removed_events)[0]

fevent=plt.figure()
# Plot most:
plt.scatter(new_robj.mw,new_robj.E_residual,s=30)
plt.scatter(new_robj.mw[event_remove_ind],new_robj.E_residual[event_remove_ind],color='red',s=30,marker='D')
plt.xlabel('Moment Magnitude')
plt.ylabel('Event term (ln residual)')
plt.title('Rejected Recordings - Event Terms')

# Save figures:
fevent.savefig(fig_dir+evname+'.png')
fevent.savefig(fig_dir+'pdfs/'+evname+'.pdf')


## STATIONS ##
# Figure basename:
stname='rejected_site'

sta_remove_ind=where(new_robj.sta == removed_stas)[0]

fsta=plt.figure()
# Plot most:
plt.scatter(new_robj.stnum,new_robj.site_terms,s=30)
plt.scatter(new_robj.stnum[sta_remove_ind],new_robj.site_terms[sta_remove_ind],color='green',s=30,marker='H')
plt.xlabel('Station Number')
plt.ylabel('Site Term (ln residual)')
plt.title('Rejected Recordings - Site Terms')

# Save figures:
fsta.savefig(fig_dir+stname+'.png')
fsta.savefig(fig_dir+'pdfs/'+stname+'.pdf')

## TOTAL ##
# Fig name:
totalname_mw='rejected_total_mw'

ftotal_mw=plt.figure()
# Plot most:
plt.scatter(new_robj.mw,new_robj.total_residual,s=30)
plt.scatter(new_robj.mw[event_remove_ind],new_robj.total_residual[event_remove_ind],color='red',s=30,marker='D')
plt.scatter(new_robj.mw[sta_remove_ind],new_robj.total_residual[sta_remove_ind],color='green',s=30,marker='H')
plt.xlabel('Moment Magnitude')
plt.ylabel('Total Residual (ln residual)')
plt.title('Rejected Recordings by Magnitude - Total Residual')

# Save figures:
ftotal_mw.savefig(fig_dir+totalname_mw+'.png')
ftotal_mw.savefig(fig_dir+'pdfs/'+totalname_mw+'.pdf')

# By distance
totalname_r='rejected_total_r'

ftotal_r=plt.figure()
# Plot most:
plt.scatter(new_robj.r,new_robj.total_residual,s=30)
plt.scatter(new_robj.r[event_remove_ind],new_robj.total_residual[event_remove_ind],color='red',s=30,marker='D')
plt.scatter(new_robj.r[sta_remove_ind],new_robj.total_residual[sta_remove_ind],color='green',s=30,marker='H')
plt.xlabel('Rrup (km)')
plt.ylabel('Total Residual (ln residual)')
plt.title('Rejected Recordings by Distance - Total Residual')

# Save figures:
ftotal_r.savefig(fig_dir+totalname_r+'.png')
ftotal_r.savefig(fig_dir+'pdfs/'+totalname_r+'.pdf')


## PATH ##
# Figure name:
pathname='rejected_path'

fpath=plt.figure()
# Plot most:
plt.scatter(new_robj.mw,new_robj.path_terms,s=30)
plt.scatter(new_robj.mw[event_remove_ind],new_robj.path_terms[event_remove_ind],color='red',s=30,marker='D')
plt.scatter(new_robj.mw[sta_remove_ind],new_robj.path_terms[sta_remove_ind],color='green',s=30,marker='H')
plt.xlabel('Moment Magnitude')
plt.ylabel('Path Terms (ln residual)')
plt.title('Rejected Recordings - Path Terms')

# Save figures:
fpath.savefig(fig_dir+pathname+'.png')
fpath.savefig(fig_dir+'pdfs/'+pathname+'.pdf')