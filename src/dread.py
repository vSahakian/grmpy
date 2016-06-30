######Data Module######
#VJS 6/2016

#Module to read and digest data of different forms, and prepare it for grmpepy

def mread(flatfile,hashfile):
    '''
    Read data from Annemarie's flatfile
    VJS 6/2016
    
    Input 
        flatfile:   String with path to the anza flatfile from AB
        hasfile:    String with path to the anza hash file from AB, with locations
    Output
    
    '''
    
    import scipy.io as sio
    import cdefs
    from numpy import genfromtxt,where,zeros
    
    
    #Read in flatfile
    datr=sio.loadmat(flatfile)
    hashr=genfromtxt(hashfile)
    
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
