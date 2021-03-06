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

what_home=1

if what_home==0:
    #Desktop:
    HOME='/media/vsahakian/katmai'
elif what_home==1:
    #Mac:
    HOME='/Users/vsahakian'

############
############

# Final database name:
dbname = 'v6anza2013_pgrid_5sta'

#Janine's flatfile:
flatfile_pga2013=HOME+'/anza/data/databases/v6anza2013/PGA2013.dat'
flatfile_pgv2013=HOME+'/anza/data/databases/v6anza2013/PGV2013.dat'

flatfile_pgaM3=HOME+'/anza/data/databases/v6anza2013/PGA_anzaM3_1998_2016.dat'
flatfile_pgvM3=HOME+'/anza/data/databases/v6anza2013/PGV_anzaM3_1998_2016.dat'

##Vs30 file:
#vs30file=HOME+'/anza/data/vs30/California.xyz'


#######  Vs30 ########
# Output station name file for this database:
vs30_sta_file = HOME + '/anza/data/vs30/stations_v6anza2013.txt'

# Vs30 proxy ID file, from R:
vs30_proxyIDfile = HOME +'/anza/data/vs30/v6anza2013_IDproxyVs30.txt'

# Vs30 conversion file:
vs30_conversionfile = HOME + '/anza/data/vs30/y16_terrain2vs30.txt'

# Vs30 proxy file, after converting ID:
vs30_proxyfile = HOME +'/anza/data/vs30/v6anza2013_proxyVs30.txt'

# Alan's vs30 csv file:
measured_vs30_csvfile = HOME + '/anza/data/vs30/yong_measuredVs30.csv'

# Output Vs30 file with measured vs30:
vs30_measuredfile = HOME + '/anza/data/vs30/v6anza2013_measuredVs30.txt'

#Output Vs30 file:
vs30name='v6anza2013_vs30.txt'
vs30_outfile=HOME+'/anza/data/vs30/'+vs30name
######################


#Save the database object paths:
dbfname_raw=HOME+'/anza/data/databases/v6anza2013/v6anza2013.pckl'
dbfname_pgrid=HOME+'/anza/data/databases/v6anza2013/v6anza2013_pgrid.pckl'
dbfname=HOME+'/anza/data/databases/v6anza2013/v6anza2013_pgrid_5sta.pckl'

#Figure paths:
fig_dir=HOME+'/anza/data/databases/v6anza2013/figs/'


# Other Parameters:
# If a station/event is outside the propagation grid, remove it:
propgrid_W = -118.12
propgrid_S = 32.38
dx = 0.01
dy = 0.01
nx = 278
ny = 218

#Number of grid points to use as buffer for keeping inside the propagation grid:
buffer=7

#Minimum number of stations for each event:
min_stations=5

print 'Will use a buffer of %f, and minimum number of stations as %f' % (buffer,min_stations)

################
##  Read data ##
################

print 'Reading in data...'

#### Read in PGA/PGV for 2013 and M greater than 3 ####

print 'PGA 2013'
# PGA 2013
evnum_pga2013,evlat_pga2013,evlon_pga2013,evdep_pga2013,sta_pga2013,stlat_pga2013,stlon_pga2013,stelv_pga2013,grcircle_pga2013,ml_pga2013,mw_pga2013,pga_millig_pga2013,source_i_pga2013,receiver_i_pga2013,pga_snr_pga2013=dr.read_jsbfile(flatfile_pga2013)

print 'PGV 2013'
# PGV 2013
evnum_pgv2013,evlat_pgv2013,evlon_pgv2013,evdep_pgv2013,sta_pgv2013,stlat_pgv2013,stlon_pgv2013,stelv_pgv2013,grcircle_pgv2013,ml_pgv2013,mw_pgv2013,pgv_cmsec_pgv2013,source_i_pgv2013,receiver_i_pgv2013,pgv_snr_pgv2013=dr.read_jsbfile(flatfile_pgv2013)

print 'PGA M3'
evnum_pgaM3,evlat_pgaM3,evlon_pgaM3,evdep_pgaM3,sta_pgaM3,stlat_pgaM3,stlon_pgaM3,stelv_pgaM3,grcircle_pgaM3,ml_pgaM3,mw_pgaM3,pga_millig_pgaM3,source_i_pgaM3,receiver_i_pgaM3,pga_snr_pgaM3,ev_year_pgaM3=dr.read_jsbfile(flatfile_pgaM3,get_year='yes')

print 'PGV M3'
evnum_pgvM3,evlat_pgvM3,evlon_pgvM3,evdep_pgvM3,sta_pgvM3,stlat_pgvM3,stlon_pgvM3,stelv_pgvM3,grcircle_pgvM3,ml_pgvM3,mw_pgvM3,pgv_cmsec_pgvM3,source_i_pgvM3,receiver_i_pgvM3,pgv_snr_pgvM3,ev_year_pgvM3=dr.read_jsbfile(flatfile_pgvM3,get_year='yes')

##
# Get the matched set of PGA for 2013:
evnum_2013,evlat_2013,evlon_2013,evdep_2013,sta_2013,stlat_2013,stlon_2013,stelv_2013,grcircle_2013,ml_2013,mw_2013,pga_millig_2013,pga_snr_2013,pgv_cmsec_2013,pgv_snr_2013,source_i_2013,receiver_i_2013 = dr.match_pga_pgv(evnum_pga2013,evlat_pga2013,evlon_pga2013,evdep_pga2013,sta_pga2013,stlat_pga2013,stlon_pga2013,stelv_pga2013,grcircle_pga2013,ml_pga2013,mw_pga2013,pga_millig_pga2013,pga_snr_pga2013,evnum_pgv2013,sta_pgv2013,pgv_cmsec_pgv2013,pgv_snr_pgv2013)

# And for the M greater than 3 dataset:
evnum_M3,evlat_M3,evlon_M3,evdep_M3,sta_M3,stlat_M3,stlon_M3,stelv_M3,grcircle_M3,ml_M3,mw_M3,pga_millig_M3,pga_snr_M3,pgv_cmsec_M3,pgv_snr_M3,source_i_M3,receiver_i_M3,evyear_M3 = dr.match_pga_pgv(evnum_pgaM3,evlat_pgaM3,evlon_pgaM3,evdep_pgaM3,sta_pgaM3,stlat_pgaM3,stlon_pgaM3,stelv_pgaM3,grcircle_pgaM3,ml_pgaM3,mw_pgaM3,pga_millig_pgaM3,pga_snr_pgaM3,evnum_pgvM3,sta_pgvM3,pgv_cmsec_pgvM3,pgv_snr_pgvM3,evyear=ev_year_pgaM3)


## Now remove events from 2013 from the M3 dataset so they're not replicated twice...
M3_not2013_indices = np.where(evyear_M3!=2013)[0]

# Now sample...
evnum_M3 = evnum_M3[M3_not2013_indices]
evlat_M3 = evlat_M3[M3_not2013_indices]
evlon_M3 = evlon_M3[M3_not2013_indices]
evdep_M3 = evdep_M3[M3_not2013_indices]
sta_M3 = sta_M3[M3_not2013_indices]
stlat_M3 = stlat_M3[M3_not2013_indices]
stlon_M3 = stlon_M3[M3_not2013_indices]
stelv_M3 = stelv_M3[M3_not2013_indices]
grcircle_M3 = grcircle_M3[M3_not2013_indices]
ml_M3 = ml_M3[M3_not2013_indices]
mw_M3 = mw_M3[M3_not2013_indices]
pga_millig_M3 = pga_millig_M3[M3_not2013_indices]
pga_snr_M3 = pga_snr_M3[M3_not2013_indices]
pgv_cmsec_M3 = pgv_cmsec_M3[M3_not2013_indices]
pgv_snr_M3 = pgv_snr_M3[M3_not2013_indices]
evyear_M3 = evyear_M3[M3_not2013_indices]

# Recompute the source and reciever indices:
source_i_M3,receiver_i_M3 = dr.find_source_receiver_indices(evnum_M3,sta_M3)


# Concatenate the years on the M3 dataset with event numbers, to be sure they're always unique:
evnum_M3_string = evnum_M3.astype('str')
evyear_M3_string = evyear_M3.astype('int').astype('str')
evnum_M3_withYear_string = np.core.defchararray.add(evyear_M3_string,evnum_M3_string)
evnum_M3_withYear = evnum_M3_withYear_string.astype('float')

# Now add these on for final arrays, to use for database:
evnum = np.r_[evnum_2013,evnum_M3_withYear]
evlat = np.r_[evlat_2013,evlat_M3]
evlon = np.r_[evlon_2013,evlon_M3]
evdep = np.r_[evdep_2013,evdep_M3]
sta = np.r_[sta_2013,sta_M3]
stlat = np.r_[stlat_2013,stlat_M3]
stlon = np.r_[stlon_2013,stlon_M3]
stelv = np.r_[stelv_2013,stelv_M3]
grcircle = np.r_[grcircle_2013,grcircle_M3]
ml = np.r_[ml_2013,ml_M3]
mw = np.r_[mw_2013,mw_M3]
pga_millig = np.r_[pga_millig_2013,pga_millig_M3]
pga_snr = np.r_[pga_snr_2013,pga_snr_M3]
pgv_cmsec = np.r_[pgv_cmsec_2013,pgv_cmsec_M3]
pgv_snr = np.r_[pgv_snr_2013,pgv_snr_M3]


# Get source_i and receiver_i:
source_i,receiver_i = dr.find_source_receiver_indices(evnum,sta)

print 'Computing Rrup'
#Compute Rrup for the data:
rrup=dr.compute_rrup(evlon,evlat,evdep,stlon,stlat,stelv)




###############
#    Vs30    ##
###############

print 'Acquiring Vs30 from grid file'
#Get Vs30:
#Find unique stations so it can be saved to a file later:
unique_stations=np.unique(sta)
sta_index=[]
for sta_i in range(len(unique_stations)):
    sta_index.append(int(np.where(sta==unique_stations[sta_i])[0][0]))

#Get the lat and lon of these points only:
stlon_unique=stlon[sta_index]
stlat_unique=stlat[sta_index]

# Write out the file to use for vs30 extraction: 
#           sta  lon lat

stationfile = open(vs30_sta_file,'w')

for station in range(len(unique_stations)):
    stationfile.write('%s\t%12.8f\t%10.8f\n' % (unique_stations[station], stlon_unique[station], stlat_unique[station]))
stationfile.close()


##Interpolate for vs30:
#vs30_unique=dr.interp_vs30(stlat_unique,stlon_unique,vs30file)
#
##Save them to a file:
#save_vs30=np.c_[unique_stations,stlon_unique,stlat_unique,vs30_unique]
#vs30fmt="%5s %9s %7s %7s"
#np.savetxt(vs30_outfile,save_vs30,fmt=vs30fmt)
#

#######
#######

###
## BEfore doing this, run /software/r/yong_vs30proxy.r, setting vs30_sta_file and vs30_proxyIDfile paths in it.
###

##  Get proxy values for every site in database:
dr.vs30proxy_id2vs30(vs30_proxyIDfile,vs30_conversionfile,vs30_proxyfile)

dr.get_measured_vs30(vs30_sta_file,measured_vs30_csvfile,vs30_measuredfile)

sta_combined,vs30_unique,vs30_method_unique = dr.combine_measured_proxy_vs30(vs30_proxyfile,vs30_measuredfile,vs30_outfile)

########
########

# But this is just for the unique stations - so what about each record?
#Now read in the vs30 for each recording:
vs30=np.zeros(len(stlon))
vs30_method=[]

for sta_i in range(len(sta)):
    #Where is the corresponding station informatioN?
    vs30ind=np.where(unique_stations==sta[sta_i])[0]
    
    vs30[sta_i]=vs30_unique[vs30ind]
    vs30_method.append(vs30_method_unique[vs30ind][0])

# Turn vs30_method into an array
vs30_method = np.array(vs30_method)

#Now vs30 cna be used to go into the database...




############################
#######  Make database  ####
############################

print 'Starting to make database...'

##Make stnum
#Make a list, of same length of the unique string station variable unique_stations,
#   with the corresponding unique numbers:
unique_stnum=np.array(range(len(unique_stations)))+1

#Initiate the stnum array, of same length as sta:
stnum=np.zeros(len(sta))

#Loop through the recordings (sta), to get the number that corresponds to them:
for sta_i in range(len(sta)):
    #Get the index of the unique stations that this recoridng corresponds to:
    unique_sta_index=np.where(unique_stations==sta[sta_i])[0]
    #Now pull the station number that corresponds to this:
    stnum[sta_i]=unique_stnum[unique_sta_index]    

print 'Converting PGA'
##   
##Convert pga_milli g to pga in m/s/s
pga=(pga_millig/1000)*9.81

##
##Set pgv - convert cm/sec to m/s
pgv=pgv_cmsec/100

print 'Make database'
##
##Make database:
db2016_1=cdf.db(evnum,sta,stnum,ml,mw,pga,pgv,rrup,vs30,evlat,evlon,evdep,stlat,stlon,stelv,source_i,receiver_i,vs30_method=vs30_method,pga_snr=pga_snr,pgv_snr=pgv_snr)

##Save database:
datobj=open(dbfname_raw,'w')
pickle.dump(db2016_1,datobj)
datobj.close()

print 'Database saved to %s' % dbfname_raw



############################
#### Print bad records #####
####     to file        ####
############################

# Some of the records have really long distances, turns out because the event
# info has latitude in both the latitude/longitude columns:
bad_evnum = evnum[np.where(rrup>500)[0]].astype('int')
bad_sta = sta[np.where(rrup>500)[0]]

# concatenate:
bad_output = np.c_[bad_evnum,bad_sta]

f = open('/Users/vsahakian/anza/data/databases/v6anza2013/badlines_pga2013.txt','w')
for linei in range(len(bad_evnum)):
    line = np.str(bad_evnum[linei]) + '\t' + np.str(bad_sta[linei]) + '\n'
    f.write(line)
f.close()



#########################
#####Sample Database#####
#########################

print 'Sampling database'

# Some recordigns are bad because the lat/lon is printed twice, so distances
#   are huge...find where this is NOT the case, and keep those indices.
#  At the same time, remove places where pga_snr or pgv_snr are equal to infinity...
rrup_snr_good_indices = np.where((rrup<500) & (pga_snr!=np.inf) & (pgv_snr!=np.inf) & (pga_snr>=4) & (pgv_snr>=2))[0]

# Sample:
dr.recording_sample(dbfname_raw,rrup_snr_good_indices,dbfname_pgrid)

##Then, sample the database so that it includes only events inside the propagation bridge
##   plus a buffer, and also recorded on a minimum number of stations:
#
#propgrid_E = propgrid_W + (dx*nx)
#propgrid_N = propgrid_S + (dy*ny)
#
## add buffer:
#propgrid = [[propgrid_W + (buffer*dx),propgrid_E - (buffer*dx)],[propgrid_S + (buffer*dy),propgrid_N - (buffer*dy)]]
#
## Sample to remove events outside of propagation grid:
#dr.db_propgrid_sample(dbfname_raw,propgrid,dbfname_pgrid)
#print 'Removed events outside of grid, saved to %s' % dbfname_pgrid

#Sample by minimum number of stations:
dr.db_station_sample(dbfname_pgrid,min_stations,dbfname)
print 'Sampled by minimum number of stations, saved to %s' % dbfname

# Set the output directory to the last one here:
dbpathout = dbfname

#######################End making database object############################


####################################################
#############Preliminary Plots######################
####################################################

print 'Making plots'

###Make plots with the sampled database###
#Open the database object;
datobj=open(dbpathout,'r')
db=pickle.load(datobj)
datobj.close()

#Plot PGA against magintude, with distance colored:
bm_min=10
bm_max=250
axlims_m=[[0.2,4.6],[-7,0]]
f_mpath=fig_dir+dbname+'_pga_mag.pdf'
f_mpng=fig_dir+dbname+'_pga_mag.png'

f_mag=db.plot_rpga(bm_min,bm_max,axlims_m)
f_mag.savefig(f_mpath)
f_mag.savefig(f_mpng)

#PGV..
axlims_m=[[0.2,4.6],[-8,-1]]
f_mpath=fig_dir+dbname+'_pgv_mag.pdf'
f_mpng=fig_dir+dbname+'_pgv_mag.png'

f_mag=db.plot_rpga(bm_min,bm_max,axlims_m,predparamflag='pgv')
f_mag.savefig(f_mpath)
f_mag.savefig(f_mpng)


######
#Plot against distance, with magnitude colored:
br_min=0.6
br_max=4.5

axlims_r=[[0,300],[-7,-1]]
f_rpath=fig_dir+dbname+'_pga_r.pdf'
f_rpng=fig_dir+dbname+'_pga_r.png'

f_r=db.plot_mpga(br_min,br_max,axlims_r)
f_r.savefig(f_rpath)
f_r.savefig(f_rpng)

# PGV:
axlims_r=[[0,300],[-8,-0]]
f_rpath=fig_dir+dbname+'_pgv_r.pdf'
f_rpng=fig_dir+dbname+'_pgv_r.png'

f_r=db.plot_mpga(br_min,br_max,axlims_r,predparamflag='pgv')
f_r.savefig(f_rpath)
f_r.savefig(f_rpng)

######
#Plot M vs. Rrup, with log10 pga colored:
bpga_min = -6.9
bpga_max = -2

axlims_pga = [[0,300],[0,5]]
f_pgapath = fig_dir+dbname+'_pga_m_rrup.pdf'
f_pgapng = fig_dir+dbname+'_pga_m_rrup.png'

f_pga = db.plot_m_rrup(bpga_min,bpga_max,axlims_pga)
f_pga.savefig(f_pgapath)
f_pga.savefig(f_pgapng)

