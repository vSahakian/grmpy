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
    return event,sta,N,ml,mw,DA,DV,dist[:,0],vs30,lat,lon,depth,stlat,stlon,source_i,receiver_i


def jb_ascii_read(flatfile):
    '''
    Read in data from Janine's test events format
    Input:
        flatfile:           String with path to the flatfile
    '''
    
    from numpy import genfromtxt,unique
    
    #First read in the stations from the flatfile:
    #Will be in two columns:   sta    chan
    star=genfromtxt(flatfile,dtype='S5',skip_header=1,usecols=[0,1])
    #Read in the data from the flatfile:
    datr=genfromtxt(flatfile,skip_header=1,usecols=range(2,12))
    
    #Find the data per station:
    #What stations are in this file?
    unique_stations=unique(star[:,0])
    
    #Loop through them, and get the horizontal components:
    





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
    
    
    
######
def db_station_sample(dbpath_in,numstas,dbpath_out):
    '''
    Sample a database to only include events recorded on a minimum number
    of stations
    VJS 8/2016
    
    Input:
        dbpath_in:          String with path to the input database
        numstas:            Minimum number of stations for events to be recorded on
        dbpath_out:         STring with path to output database
    Output:
        Writes out the sampled database to dbpath_out
    '''
    import cPickle as pickle
    from numpy import unique,where,array,r_
    import cdefs as cdf
    
    #Read in the original database
    dbfile=open(dbpath_in,'r')
    db_orig=pickle.load(dbfile)
    dbfile.close()
    
    #Find how many unique events there are:
    unique_events=unique(db_orig.evnum)
    nevents=len(unique_events)
    
    #Initiate the "keep" event index array:
    keep_event_ind=array([]).astype('int')
    
    #Loop over the unique events and determine if they are recorded on the 
    #minimum number of stations:
    
    for event_ind in range(nevents):
        #Call the event:
        eventi=unique_events[event_ind]
        #Find where in the database there are recordings of this event:
        db_event_ind=where(db_orig.evnum==eventi)[0].astype('int')
        
        #Get the stations for this event:
        eventi_stas=db_orig.sta[db_event_ind]
        #Get the unique stations recording this event:
        num_unique_stas_i=len(eventi_stas)
        
        #If it's at least the number of minimum stations, keep this stuff, save 
        #it to the keep index variable:
        if num_unique_stas_i>=numstas:
            keep_event_ind=r_[keep_event_ind,db_event_ind]
    
    
    #Save the keep event variable as integers:
    keep_event_ind.astype('int')
    
    #Now save just these indices in the database:
    edepth=db_orig.edepth[keep_event_ind]
    elat=db_orig.elat[keep_event_ind]
    elon=db_orig.elon[keep_event_ind]
    evnum=db_orig.evnum[keep_event_ind]
    ffdf=db_orig.ffdf[keep_event_ind]
    md_ffdf=db_orig.md_ffdf[keep_event_ind]
    ml=db_orig.ml[keep_event_ind]
    mw=db_orig.mw[keep_event_ind]
    pga=db_orig.pga[keep_event_ind]
    pga_pg=db_orig.pga_pg[keep_event_ind]
    pgv=db_orig.pgv[keep_event_ind]
    r=db_orig.r[keep_event_ind]
    receiver_i=db_orig.receiver_i[keep_event_ind]
    source_i=db_orig.source_i[keep_event_ind]
    sta=db_orig.sta[keep_event_ind]
    stlat=db_orig.stlat[keep_event_ind]
    stlon=db_orig.stlon[keep_event_ind]
    stnum=db_orig.stnum[keep_event_ind]
    vs30=db_orig.vs30[keep_event_ind]
    
    #BEFORE SAVING:
    #cdefs only takes DA and DV in nm/s/s and nm/s...convert to these (currently
    #in m/s/s and m/s)
    DA=pga/1e-9
    DV=pga/1e-9
    
    #Make sampled database:
    db_samp=cdf.db(evnum,sta,stnum,ml,mw,DA,DV,r,vs30,elat,elon,edepth,stlat,stlon,source_i,receiver_i)
    
    #Save to file...
    doutfile=open(dbpath_out,'w')
    pickle.dump(db_samp,doutfile)
    doutfile.close()
    