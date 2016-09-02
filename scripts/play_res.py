import dread as dr
import cdefs as cdf
import numpy as np
import inversion as inv
import pickle
import matplotlib.pyplot as plt
import gmpe as gm
import rescomp as rcomp


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

#Location to store figures:
fig_dir=HOME+'/anza/models/figs/'
obj_dir=HOME+'/anza/models/pckl/'

######Making the database object - once made and saved, do not run this again######
#
##AB's flatfile:
#ffile=HOME+'/anza/data/Anzadata_Acc_Vel_May2016_40_50.mat'
#hfile=HOME+'/anza/data/YSH_2010.hash'
#
##Read data from .mat file, store important values:
#ev,sta,N,ml,mw,DA,DV,r,vs30,lat,lon,depth=dr.mread(ffile,hfile)
#
##Load into database object
#abdb=cdf.db(ev,sta,N,ml,mw,DA,DV,r,vs30,lat,lon,depth)
#
##Save the database object:
#fname=HOME+'/anza/data/abdb.pckl'
#datobj=open(fname,'w')
#pickle.dump(abdb,datobj)
#datobj.close()
########################End making database object############################


##########Open the database object:###################
#Filename:
fname=HOME+'/anza/data/abdb.pckl'
datobj=open(fname,'r')
abdb=pickle.load(datobj)
datobj.close()

##Residual Computation####
#Load in the model to use:
fname=obj_dir+'regr_0.0_6.5_resid_2676.06963031.pckl'
model=pickle.load(open(fname,'r'))

#Overall residual, 
#In some places, vs30 is 0.  Set these to vref.
vref=760
mdep_ffdf=0
#Where are they 0?
vs30_0ind=np.where(abdb.vs30==0)[0]
#Get vs30 from database...
vs30=abdb.vs30
#Set 0 entries to vref:
vs30[vs30_0ind]=vref


#####Get Total Residuals######

#Now compute the predicted value of PGA...
d_predicted=gm.compute_model(model.m,model.rng,abdb.mw,abdb.r,abdb.ffdf,vs30,vref,mdep_ffdf)

#Get residuals:
total_residuals,mean_residual,std_dev=rcomp.total_residual(abdb,d_predicted)


###########

###Get Event and Within-Event Residuals##

##Get unique events:
unique_events=np.unique(abdb.evnum)

##Make a class for each event; append them to a list, made empty here:
event_list=[]
d_predicted_list=[]

##Zero out arrays for info used for plotting, for all events...
E_evnum=[]
E_mw=[]
E_residual=[]
E_std_dev=[]

##Zero out lists to use for indexing - use this resulting index, the station index,
# to sort out the events for the station object later on...
unique_stnums=np.unique(abdb.stnum)
unique_stas=np.unique(abdb.sta)

#Loop through the unique events, make each into an object, append to event list
for i in range(len(unique_events)):
    unique_ind=np.where(unique_events[i]==abdb.evnum)[0]
    #Get the predicted data for this event only, for each recording
    #d_predicted_i is an array, with len = # of recordings for event):
    d_predicted_i=d_predicted[unique_ind]
    
    #Get the database info for this event:
    evnum_i=abdb.evnum[unique_ind]
    sta_i=abdb.sta[unique_ind]
    stnum_i=abdb.stnum[unique_ind]
    ml_i=abdb.ml[unique_ind]
    mw_i=abdb.mw[unique_ind]
    pga_i=abdb.pga[unique_ind]
    pgv_i=abdb.pgv[unique_ind]
    pga_pg_i=abdb.pga_pg[unique_ind]
    r_i=abdb.r[unique_ind]
    ffdf_i=abdb.ffdf[unique_ind]
    md_ffdf_i=abdb.md_ffdf[unique_ind]
    lat_i=abdb.lat[unique_ind]
    lon_i=abdb.lon[unique_ind]
    depth_i=abdb.depth[unique_ind]
    
    #for sta_ind in len(unique_stnums):
    #    
    
    #The only one that is not directly from teh database is vs30; this is because
    #there were some 0's in it, so above that is compensated for by changing 
    #them to vref.
    vs30_i=vs30[unique_ind]
    
    
    #Make the event object:
    eventi=cdf.event(evnum_i,sta_i,stnum_i,ml_i,mw_i,pga_i,pgv_i,pga_pg_i,r_i,vs30_i,ffdf_i,md_ffdf_i,lat_i,lon_i,depth_i)
    
    #Compute the event terms:
    evnum_i,evmw_i,E_residual_i,std_dev_i=rcomp.event_residual(eventi,d_predicted_i)
    
    #Add the residual information to the event object:
    eventi.add_E_resid(E_residual_i,std_dev_i)
    
    #Get the Within-Event Residuals:
    evnum_i,evmw_i,sta_i,stnum_i,W_residuals_i,W_mean_i,W_std_dev_i=rcomp.within_event_residual(eventi,d_predicted_i,eventi.E_residual)    
    
    #Add the within-event residuals to the event object:
    eventi.add_W_resids(W_residuals_i,W_mean_i,W_std_dev_i)
    
    #Append the event object, and the d_predicted to the list:
    event_list.append(eventi)
    d_predicted_list.append(d_predicted_i)
    
    #Append to the residual arrays, to use later for plotting:
    E_evnum.append(evnum_i)
    E_mw.append(evmw_i)
    E_residual.append(E_residual_i)
    
#Turn those into arrays:
E_evnum=np.array(E_evnum)
E_mw=np.array(E_mw)
E_residual=np.array(E_residual)
E_std_dev=np.array(E_std_dev)
 
 
 
###Make Station Objects, for other residuals####

#Start an empty list, store all station objects here in the end:    
station_list=[]

#First loop over the list of unique stations    
for sta_ind in range(len(unique_stnums)):
    #Station number?
    station_num_i=unique_stnums[sta_ind]
    
    #Station name:
    station_name_i=unique_stas[sta_ind]
    
    #Zero out the lists/arrays that will be used to add to the station object for htis statioN:
    evnum=[]
    ml=[]
    mw=[]
    pga_pg=[]
    pga=[]
    pgv=[]
    r=[]
    ffdf=[]
    md_ffdf=[]
    lat=[]
    lon=[]
    depth=[]
    #Event residuals:
    E_residual=[]
    #Within-event residuals:
    W_residual=[]
    
    #Does this e
    for event_ind in range(len(unique_events)):
        #What event is this?
        eventi=event_list[event_ind]
        #What stations record thsi event?
        event_sta_i=eventi.stnum
        #Is the station being referenced int eh outer loop contained here?
        sta_ev_ind=np.where(event_sta_i==station_num_i)[0]
        
        #If this station records thsi event, grab the information:
        if sta_ev_ind.size!=0:
            print sta_ev_ind.size
            vs30=eventi.vs30[sta_ev_ind]
            evnum.append(eventi.evnum[sta_ev_ind])
            ml.append(eventi.ml[sta_ev_ind])
            mw.append(eventi.mw[sta_ev_ind])
            pga_pg.append(eventi.pga_pg[sta_ev_ind])
            pgv.append(eventi.pgv[sta_ev_ind])    
            r.append(eventi.r[sta_ev_ind])    
            ffdf.append(eventi.ffdf[sta_ev_ind])
            md_ffdf.append(eventi.md_ffdf[sta_ev_ind])
            lat.append(eventi.lat[sta_ev_ind])
            lon.append(eventi.lon[sta_ev_ind])
            depth.append(eventi.depth[sta_ev_ind])
            E_residual.append(eventi.E_residual)
            W_residual.append(eventi.W_residuals[sta_ev_ind])
        elif sta_ev_ind.size==0:
            continue
            
        
    #AFter looping over all events, convert these lists into arrays, to put
    #into a station object:
    #Vs30 stays as just one number
    evnum=np.array(evnum)
    ml=np.array(ml)
    mw=np.array(mw)
    pga_pg=np.array(pga_pg)
    pgv=np.array(pgv)
    r=np.array(r)
    ffdf=np.array(ffdf)
    md_ffdf=np.array(md_ffdf)
    lat=np.array(lat)
    lon=np.array(lon)
    depth=np.array(depth)
    E_residual=np.array(E_residual)
    W_residual=np.array(W_residual)
    
    #Put into a station object...
    station_i=cdf.station(station_name_i[0],station_num_i,vs30[0],evnum,ml,mw,pga_pg,pga,pgv,ffdf,md_ffdf,lat,lon,depth,E_residual,W_residual)
    
    #Append to the station list...
    station_list.append(station_i)


###Plotting####
#Plot the event residuals for each station in a different color...

#Axes:
axxlim=[1,4]
axylim=[-4,4]

#Get the range for the colorbar:
crangemax=len(station_list)
crange=np.array(range(len(station_list))).astype(float)/crangemax

#Get the colorbar:
colors=plt.cm.rainbow(crange)

plt.figure()
for station_ind in range(len(station_list)):
    #Get the station 
    station_i=station_list[station_ind]
    
    #Get what you're plotting...
    mw=station_i.mw
    W_residuals=station_i.W_residual
    
    #Get the color info and label info:
    color_i=colors[station_ind]
    sta_lab=station_i.sta
    
    print color_i
    
    #Plot
    plt.scatter(mw,W_residuals,edgecolors=color_i,facecolors='none',lw=0.8,label=sta_lab)
    
#Add legend:
plt.legend(loc=4)

#Add titles, limits...
plt.xlim(axxlim)
plt.ylim(axylim)
plt.xlabel('Moment Magnitude')
plt.ylabel('ln Residual')
plt.title('Within-Event Residuals by station')

#Show...
plt.show()




