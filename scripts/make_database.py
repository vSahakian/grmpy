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
    
    
    
#####Making the database object - once made and saved, do not run this again######

#AB's flatfile and hasfhile; Debi's station file; columns for Debi's stationfile:
ffile=HOME+'/anza/data/Anzadata_Acc_Vel_May2016_40_50.mat'
hfile=HOME+'/anza/data/YSH_2010.hash'
sfile=HOME+'/anza/data/stations/AZ_stations_elevation.ll'
scols=[0,3,4,5]

#Read data from .mat file, store important values:
ev,sta,N,ml,mw,DA,DV,r,vs30,lat,lon,depth,stlat,stlon,stelv,source_i,receiver_i=dr.mread(ffile,hfile,sfile,scols)

#Convert DA and DV to m/s/s and m/s, from nm/s/s and nm/s:
DAm=DA*1e-9
DVm=DV*1e-9

#Load into database object
abdb=cdf.db(ev,sta,N,ml,mw,DAm,DVm,r,vs30,lat,lon,depth,stlat,stlon,stelv,source_i,receiver_i)

#Save the database object:
fname=HOME+'/anza/data/abdb.pckl'
datobj=open(fname,'w')
pickle.dump(abdb,datobj)
datobj.close()

#Then, sample the database so that it includes only events recorded on a minimum
#number of stations:
min_stations=5
dbpathin=HOME+'/anza/data/abdb.pckl'
dbpathout=HOME+'/anza/data/abdb_5sta.pckl'

#Sample:
dr.db_station_sample(dbpathin,min_stations,dbpathout)

#######################End making database object############################