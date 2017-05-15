#############################
###  Database statistics  ###
#############################
##
## Module to compute 
## statistics for databases
## VJS 5/2017
##

def compute_lagdistances(sta,stnum,lon,lat):
    '''
    Compute the lag distances between all stations in the given set
    Input:
        sta:            Array with strings of station names (n x 1)
        stnum:          Array with unique station numbers (n x 1)
        lon:            Array with station longitudes (n x 1)
        lat:            Array with station latitudes (n x 1)
    Output:
        lagdistance:   Upper triangular matrix with lag distances for all station pairs (n x n)
    '''
    
    from pyproj import Geod
    import numpy as np
    
    # Make projection:
    p = Geod(ellps='WGS84')
    
    # Make lag matrix:
    lagdistance_full = np.zeros((len(stnum),len(stnum)))
    
    ## Start to fill in lag matrix
    # Loop over all stations, make a matrix with the lon and lat of just that station, and compute the distance to all other stations (lon,lat):
    for stationi in range(len(stnum)):
        azimuthi,backazimuthi,distancei = p.inv(lon[stationi]*np.ones(len(stnum)),lat[stationi]*np.ones(len(stnum)), lon, lat)
        
        # Fill the matrix with these distances for this station:
        lagdistance_full[stationi,:] = distancei/1000
    
    # Turn it into an upper triangular matrix:
    lagdistance = np.triu(lagdistance_full)    
                
    # Return lag distance:
    return lagdistance
        

##################################
##################################

def semivar_rbin(robj,bins):
    '''
    Sort the residuals object into bins, and save what's necessary for semivariogram.  
    For each variable, outputs a list of arrays, where each array is the 
        number of recordings in that bin, and the list length is the number of
        bins.  Not all bins have the same number of recordings in them.
    Input:
        robj:           Residuals object 
        bins:           Array of bin edges, i.e. [0,3,6,7,..].  number of bins = n = len(bins) - 1
    Output:
        evnum:          List of arrays with event numbers of all recordings in this bin n x (m x 1)
        sta:            List of arrays with all stations in this bin   n x (m x 1)
        stnum:          List of arrays with all station numbers in this bin, rearranged so len(unique(stnum)) = len(unique(sta)), not len(unique(sta_original) n x (m x 1)
        stlon:          List of arrays with station longitude for bin   n x (m x 1)
        stlat:          List of arrays with station latitude for bin   n x (m x 1)
        path_terms:     List of arrays with path terms for bin   n x (m x 1)
        rrup:           List of arrays with rrup for bin   n x (m x 1)
    '''
    
    import numpy as np
    
    # First, digitize the residuals object:
    digitize_data = np.digitize(robj.r,bins)
    
    # Inititate final lists to be empty, then append tot hem later:
    evnum = []
    sta = []
    stnum = []
    stlon = []
    stlat = []
    path_terms = []
    rrup = []
    
    # Loop over bins to get the necessary arrays:
    for bini in range(len(bins)-1):
        print 'Bin %i: bin[%.1f] <= x < bin[%.1f]' % (bini,bins[bini],bins[bini+1])
    
        # Get the raw output information for this bin - get stnum last:
        # Indices for this bin:
        digitize_ind = np.where(digitize_data==bini+1)[0]
           
        # If there are recordings in this bin:
        if len(digitize_ind)!=0:
            print '\t %i recordings in this bin' % len(digitize_ind)
            
            # Data for this bin:
            evnumi = robj.evnum[digitize_ind]
            stai = robj.sta[digitize_ind]
            stloni = robj.stlon[digitize_ind]
            stlati = robj.stlat[digitize_ind]
            path_termsi = np.array(robj.path_terms)[digitize_ind]
            rrupi = robj.r[digitize_ind]
        
            # stnum needs to start at 0, and be the length of the unique stations in this bin:
            unique_sta,usta_ind = np.unique(stai,return_index=True)
            # unique ones start at 0...
            unique_stnum = np.arange(0,len(unique_sta),1)
            
            # Now make the full array, and match it up to the station names:
            stnumi = np.zeros(len(stai))
            
            for stationi in range(len(stai)):
                ustation_ind = np.where(unique_sta==stai[stationi])[0]
                
                stnumi[stationi] = unique_stnum[ustation_ind]
    
        elif len(digitize_ind)==0:
            print '\t No recordings in this bin, adding empty arrays for it'
            evnumi = np.array([])
            stai = np.array([])
            stnumi = np.array([])
            stloni = np.array([])
            stlati = np.array([])
            path_termsi = np.array([])
            rrupi = np.array([])
            
        #Now concatenate them for the final output:
        evnum.append(evnumi)
        sta.append(stai)
        stnum.append(stnumi)
        stlon.append(stloni)
        stlat.append(stlati)
        path_terms.append(path_termsi)
        rrup.append(rrupi)

    # Return output:
    return evnum,sta,stnum,stlon,stlat,path_terms,rrup
    
    
##################################
##################################
def semivar_variables(evnum,sta,stnum,stlon,stlat,path_terms,rrup,lagdistance_x):
    '''
    Compute and return the lag distance (x) and semivariance (y) variables 
    for the semivariogram plot.
    Input:
        evnum:              Array with event numbers per recording in semivariogram bin (m x 1)
        sta:                Array with station codes per recording in semivariogram bin (m x 1)
        stnum:              Array with station numbers per recording in semivariogram bin (m x 1)
        stlon:              Array with station longitudes per recording in semivariogram bin (m x 1)
        stlat:              Array with station latitudes per recording in semivariogram bin (m x 1)
        path_terms:         Array with path terms per recording in semivariogram bin (m x 1)
        rrup:               Array with rrup per recording in semivariogram bin (m x 1)
        lagdistance_x:      Array with lag distance; also defines the bins such that bin[n] <= station_separation < bin[n+1]
    Output:
        lagdistance_x:      Same as input - Array with lag distance variable for this semivariogram bin
        semivariance_y:     Array with semivariance for this bin
    '''

    import numpy as np
    
    
    # Compute station to station distance for all stations involved in this bin:
    print 'Computing station to station distance for all stations involved in this semivariogram'
    
    # get unique stations...
    unique_sta,ustaind = np.unique(sta,return_index=True)
    unique_stnum = stnum[ustaind]
    unique_stlon = stlon[ustaind]
    unique_stlat = stlat[ustaind]
    
    # get lag distance matrix:
    lagdistance_matrix = compute_lagdistances(unique_sta,unique_stnum,unique_stlon,unique_stlat)
    
    # Get unique earthquakes:
    unique_evnum,uevind = np.unique(evnum,return_index=True)
    
    # Initiate the final semivariance variable:
    semivariance_y = np.array([])
    
    print 'Looping over lag distance bins and earthquakes'
    # Loop over lag distance bins:
    for lagbin_i in range(len(lagdistance_x)-1):
        # Bin min and max:
        lagbin_min_i = lagdistance_x[lagbin_i]
        lagbin_max_i = lagdistance_x[lagbin_i+1]
        
        # If this is the first bin, and the minimum (left) bin value is 0, make it greater than 0, but small:
        if lagbin_min_i==0:
            lagbin_min_i = 1/10000000.
            
        # For this lag bin, initiate an empty array for the square of the path term differences:
        path_term_diff_square_bin = np.array([])
        
        # Now, loop over every earthquake...
        for earthquake_i in range(len(unique_evnum)):
            
            # Which stations are recorded by this earthquake?
            evnum_ind = np.where(evnum==unique_evnum[earthquake_i])[0]
            stnum_list_eqi = stnum[evnum_ind].astype('int')
            
            # Get the path terms for just this earthquake:
            path_terms_eqi = path_terms[evnum_ind]
            
            # Use stnum_list_eqi as the Rows and Columns array to make new station2statoin lag distnace matrix:
            # first keep the columns -
            lagdistance_columns = lagdistance_matrix[:,stnum_list_eqi]
            # then chop to keep the rows, and get the final sampled lag distance matrix:            
            lagdistance_matrix_eqi = lagdistance_columns[stnum_list_eqi,:]
            
        
            # Find which stations recording this earthquake are in this lag bin, if any:
            stations_in_lagbin_ind = np.where((lagdistance_matrix_eqi >= lagbin_min_i) & (lagdistance_matrix_eqi < lagbin_max_i)) 
            # rows:
            r_i = stations_in_lagbin_ind[0]
            c_i = stations_in_lagbin_ind[1]
            
            # If there are at least two stations (two paths) in this bin (one distance),or more, compute path term differneces:
            if len(r_i)>=1:
                path_term_diff_square_eqi =  (path_terms_eqi[r_i] - path_terms_eqi[c_i])**2
            
                # Append to the list of path term difference squares for this lag bin:
                path_term_diff_square_bin = np.r_[path_term_diff_square_bin,path_term_diff_square_eqi]
            else:
                continue
            
        # After looping over all earthquakes, but still in the lag bin, sum the squares of differneces, and divide by 2*N
        # - but only if there are more than one entry in this bin:
        if len(path_term_diff_square_bin) >=1:
            semivariance_bin_i = np.sum(path_term_diff_square_bin)/(2*len(path_term_diff_square_bin))
        
        # Otherwise, set the value to 'NaN':
        else:
            semivariance_bin_i = np.nan
        
        # Concatenate to the end of this bin:
        semivariance_y = np.r_[semivariance_y,semivariance_bin_i]

                
    # Return values after looping over all bins:
    return lagdistance_x, semivariance_y
     

##################################
##################################
def subplot_all_semivariograms(lagdistance_x_all,semivariance_y_all,ncols,axlimits,semivar_bins):
    '''
    Plot all semivariograms at once on a subplot
    Input:
        lagdistance_x_all:          List of arrays of lag distance (x) to plot
        semivariance_y_all:         List of arrays of semivariance (y) to plot
        ncols:                      Number of columns for subplot
        axlimits:                   Axis limits for plot: [[xmin,xmax],[ymin,ymax]]
        semivar_bins:               Array with semivariance rrup bins
    Output:
        semivariograms_fig:         Figure with semivariogram subplots
    '''
    
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator
    from numpy import ceil,remainder,where,isnan,floor
    
    # How many semivariograms/rrup bins?
    nbins = len(semivariance_y_all)
    
    # How many rows?
    nrows = ceil(nbins/ncols)
    
    # Initiate figure:
    semivariograms_fig, axarr = plt.subplots(int(nrows),int(ncols),figsize=(5*ncols,3*nrows))
    
    # Loop over bins:
    for semivar_bin in range(nbins):
        # Get row number:
        row_i = floor(semivar_bin/ncols)
        # Get column number:
        col_i = remainder(semivar_bin,ncols)
        
        print 'row %i, column %i' % (row_i,col_i)
        
        ax = axarr[row_i,col_i]
        
        x = lagdistance_x_all[semivar_bin]
        y = semivariance_y_all[semivar_bin]
        
        # Find where it's not a nan:
        keep_i = where(isnan(y)==False)[0]
        x = x[keep_i]
        y = y[keep_i]
        
        # Scatter:
        ax.scatter(x,y,facecolors='#69bcb7',edgecolors='k',s=25,label='Rrup bin %.1f to %.1f km' % (semivar_bins[semivar_bin],semivar_bins[semivar_bin+1]))
        
        # Axis limits:
        ax.set_xlim(axlimits[0])
        ax.set_ylim(axlimits[1])
        
        # Set y label spacing:
        ymajorLocator = MultipleLocator(0.1)
        ax.yaxis.set_major_locator(ymajorLocator)
        
        # Set y tick spacing
        yminorLocator = MultipleLocator(0.02)
        ax.yaxis.set_minor_locator(yminorLocator)
        
        # Set y tick length:
        ax.tick_params(which='major',length=7,width=1)
        ax.tick_params(which='minor',length=4,width=1)
        
        # Labels:
        ax.set_xlabel('Lag Distance (km)')
        ax.set_ylabel('Path Semivariance')
        
        # Put the title inside:
        #ax.annotate('Rrup bin %.1f km to %.1f km' % (semivar_bins[semivar_bin],semivar_bins[semivar_bin+1]), fontsize=12, xy=())
        
        # Title/legend:
        ax.legend(frameon=False,fontsize=12,markerscale=0,loc=2)
        
    # Adust width between plots so title fits:
    plt.subplots_adjust(left=0.1,right=0.95,bottom=0.05,top=0.95,wspace=0.3,hspace=0.4)
   
    return semivariograms_fig
        