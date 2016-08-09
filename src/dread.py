######Data Module######
#VJS 6/2016

#Module to read and digest data of different forms, and prepare it for grmpepy

def mread(flatfile,hashfile,stationfile):
    '''
    Read data from Annemarie's flatfile
    VJS 6/2016
    
    Input 
        flatfile:        String with path to the anza flatfile from AB
        hashfile:        String with path to the anza hash file from AB, with locations
        stationfile:     STring with path to the station file, for anza stations 
    Output
    
    '''
    
    import scipy.io as sio
    import cdefs
    from numpy import genfromtxt,where,zeros,unique
    
    
    #Read in flatfile
    datr=sio.loadmat(flatfile)
    hashr=genfromtxt(hashfile)
    #Read in station info from columns of station file:
    #Station name:
    stat_sta_r=genfromtxt(stationfile,skip_header=3,usecols=1,dtype='S5')
    stat_coord_r=genfromtxt(stationfile,skip_header=3,usecols=[2,3])
    
    #Info from the hashfile:
    hevent=hashr[:,6]
    hlat=hashr[:,7]
    hlon=hashr[:,8]
    hdepth=hashr[:,9]
    
    #Extract components:
    devent=datr['event']            #event number
    dsta=datr['sta']                #station name
    dhdrs=datr['hdrs']              #header info...see below...
    dN=datr['N']                    #station number         
    dMl=datr['Ml']                  #catalog/local mag.
    dMw=datr['Mw']                  #Ml converted to Mw
    dDA=datr['DA']                  #PGA, geometrically averaged from 2 horiz. comp.
    dDV=datr['DV']                  #PGV, geometrically averaged from 2 horiz. comp.
    dPGA=datr['PGA']                #DA, adjusted for Vs30 and converted to g
    dPGV=datr['PGV']                #DV, adjusted for Vs30 and converted to g
    dlogPGA=datr['logPGA']          #DA adjusted to g
    dlogPGV=datr['logPGV']          #DV adjusted to cm/s
    dnewlogPGA=datr['newlogPGA']    #log10 of PGA adjusted to 10km
    dnewlogPGV=datr['newlogPGV']    #log10 of PGV adjusted to 10km
    dVs30=datr['Vs30']              #Vs30 for each station, reference in paper.
    
    #Find the lat/lon that corresponds to this event:
    #Zero out lat/lon arrays:
    ev_lat=zeros((len(devent[0]),1))
    ev_lon=zeros((len(devent[0]),1))
    ev_dep=zeros((len(devent[0]),1))
    
    #Find which row (called event_ind) in the hashfile corresponds to this event:
    for i in range(len(devent[0])):
        event_ind=where(hevent==devent[0,i])[0][0]
        #Set that row in lat/lon to this value:
        ev_lat[i]=hlat[event_ind]
        ev_lon[i]=hlon[event_ind]
        ev_dep[i]=hdepth[event_ind]
    
    #Find row from station data that corresponds to the station:
    #First zero out arrays:
    sta_lat=zeros((len(dsta[0]),1))
    sta_lon=zeros((len(dsta[0]),1))
    
    for j in range(len(dsta[0])):
        stationi=dsta[0][j][0].astype('S5')
        #Find which index of the station file this corresponds to for this station:
        station_ind=where(stationi==stat_sta_r)[0]
        
        #Now get the information for this station, lat and lon:
        sta_lat[j]=stat_coord_r[station_ind,0]
        sta_lon[j]=stat_coord_r[station_ind,1]
        
    
    ###Get indices for event and station for the sources.in and receivers.in files####
    ##Events first:
    
    #Get the unique station and event indices:
    unique_events=unique(devent)
    
    #Zero out source ind array:
    source_ind=zeros((len(devent[0])))
    #For each event in the record, devent, give it the source index to be used:
    for event_ind in range(len(unique_events)):
        eventi=unique_events[event_ind]
        
        #Find where in the recordings list the event number is the same as this one:
        recording_event_ind=where(devent[0]==eventi)
        
        #Set the source ind to be one plus this event, so it indexes with the raytracing program:
        source_ind[recording_event_ind]=event_ind+1
    
    #Now set these to integers...
    source_ind=source_ind.astype('int64')
    
    ##Next stations:
    unique_stations=unique(dsta)
    
    #Zero out array:
    receiver_ind=zeros((len(dsta[0])))
    #Loop through the unique stations:
    for station_ind in range(len(unique_stations)):
        stationi=unique_stations[station_ind]
        
        #Find where in the recordings list the station is the same as this one:
        recording_station_ind=where(dsta[0]==stationi)[0]
    
        #Set the receiver ind to be one plus this station, so it indexes with the raytracin gprogram:
        receiver_ind[recording_station_ind]=station_ind+1    
        
    #Set these to integers:
    receiver_ind=receiver_ind.astype('int64')
  
        
    #Put into arrays (instead of array of arrays):
    event=devent[0]
    sta=dsta[0]
    N=dN[0]
    ml=dMl[:,0]
    mw=dMw[:,0]
    DA=dDA[:,0]
    DV=dDV[:,0]
    PGA=dPGA[:,0]
    PGV=dPGV[:,0]
    logPGA=dlogPGA[:,0]
    logPGV=dlogPGV[:,0]
    newlogPGA=dnewlogPGA[:,0]
    newlogPGV=dnewlogPGV[:,0]
    vs30=dVs30[0]
    lat=ev_lat[:,0]
    lon=ev_lon[:,0]
    depth=ev_dep[:,0]
    stlat=sta_lat[:,0]
    stlon=sta_lon[:,0]
    source_i=source_ind
    receiver_i=receiver_ind
    
    ##
    #Split apart the header (dhdrs)
    ddepmin=dhdrs['depmin']         #dependent (y) min
    ddepmax=dhdrs['depmax']         #dependent (y) max
    dmag=dhdrs['mag']               #magnitude (catalog)
    ddist=dhdrs['dist']             #distance from source to site, epicentral...
    daz=dhdrs['az']                 #source-site azimuth
    dyear=dhdrs['year']             #year
    dday=dhdrs['day']               #day
    dhour=dhdrs['hour']             #hour
    dmin=dhdrs['min']               #min
    dsec=dhdrs['sec']               #sec
    dmsec=dhdrs['msec']             #msec
    dnevid=dhdrs['nevid']           #event number
    didep=dhdrs['idep']             #units of the independent variable:
                                    #5 implies no instrument response removed
                                    #7 is velocity in nm/s
                                    #8 is acceleration in nm/s/s
    
    #Put into useful values...
    depmin=ddepmin[0][0]
    depmax=ddepmax[0][0]
    mag=dmag[0][0]
    dist=ddist[0][0]
    az=daz[0][0]
    year=dyear[0][0]
    day=dday[0][0]
    minu=dmin[0][0]
    sec=dsec[0][0]
    msec=dmsec[0][0]
    nevid=dnevid[0][0]
    idep=didep[0][0]
    
    #Return the event numeber, station name, station number, ml,mw, PGA,pgv, 
    #epcentral distance (Dist), vs30
    return event,sta,N,ml,mw,DA,DV,dist[:,0],vs30,lat,lon,depth


def read_obj_list(obj):
    '''
    Read a pickle file with a list of event of station objects
    Input:
        obj:        String with path to the pickle file containing the list of objects
    Output:
        obj_list:   List with event or station objects in it
    '''
    
    import cPickle as pickle
    
    #Open the file
    obj=open(obj,'r')
    
    #Zero out a list to add the event objects to:
    obj_list=[]
    
    #REad the data...
    readfile=True
    while readfile==True:
        try:
            obji=pickle.load(obj)
            obj_list.append(obji)
        except EOFError:
            print ('File read complete')
            readfile=False
            obj.close()
        
    return obj_list
    