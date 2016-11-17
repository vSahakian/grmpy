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
flatfile=HOME+'/anza/data/databases/db2013_test/PGA_2013.dat'
#Vs30 file:
vs30file=HOME+'/anza/data/vs30/California.xyz'

#Output Vs30 file:
vs30name='vs30_2013test.txt'
vs30_outfile=HOME+'/anza/data/vs30/'+vs30name

#Save the database object paths:
dbfname_raw=HOME+'/anza/data/databases/db2013_test/db2013test.pckl'
dbfname=HOME+'/anza/data/databases/db2013_test/db2013test_5sta.pckl'

#Figure paths:
fig_dir=HOME+'/anza/data/databases/db2013_test/figs/'


#############
#############

print 'Reading in data'
#Read in the data:
evnum,evlat,evlon,evdep,sta,stlat,stlon,stelv,grcircle,ml,mw,pga_millig,source_i,receiver_i=dr.read_jsbfile(flatfile)

print 'Computing Rrup'
#Compute Rrup for the data:
rrup=dr.compute_rrup(evlon,evlat,evdep,stlon,stlat,stelv)

######################
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

#Interpolate for vs30:
vs30_unique=dr.interp_vs30(stlat_unique,stlon_unique,vs30file)

#Save them to a file:
save_vs30=np.c_[unique_stations,stlon_unique,stlat_unique,vs30_unique]
vs30fmt="%5s %9s %7s %7s"
np.savetxt(vs30_outfile,save_vs30,fmt=vs30fmt)

#Now read in the vs30 for each recording:
vs30=np.zeros(len(stlon))
for sta_i in range(len(sta)):
    #Where is the corresponding station informatioN?
    vs30ind=np.where(unique_stations==sta[sta_i])[0]
    vs30[sta_i]=vs30_unique[vs30ind]

#Now vs30 cna be used to go into the database...

#######################
######Make database####
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
##Set pgv to zeros, since it doesn't exist for now...
pgv=np.zeros(len(pga))

print 'Make database'
##
##Make database:
db2013test=cdf.db(evnum,sta,stnum,ml,mw,pga,pgv,rrup,vs30,evlat,evlon,evdep,stlat,stlon,stelv,source_i,receiver_i)

##Save database:
datobj=open(dbfname_raw,'w')
pickle.dump(db2013test,datobj)
datobj.close()

print 'Database saved'

########################
####Sample Database#####
########################

print 'Sampling database'
#Then, sample the database so that it includes only events recorded on a minimum
#number of stations:
min_stations=5
dbpathin=dbfname_raw
dbpathout=dbfname

#Sample:
dr.db_station_sample(dbpathin,min_stations,dbpathout)

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

#Plot against magintude, with distance colored:
bm_min=10
bm_max=140
axlims_m=[[0,4.5],[-8,-1]]
f_mpath=fig_dir+'db2013test_mag.pdf'

f_mag=db.plot_rpga(bm_min,bm_max,axlims_m)
f_mag.savefig(f_mpath)

#Plot against distance, with magnitude colored:
br_min=0
br_max=3

axlims_r=[[0,260],[-8,-1]]
f_rpath=fig_dir+'db2013test_r.pdf'

f_r=db.plot_mpga(br_min,br_max,axlims_r)
f_r.savefig(f_rpath)
