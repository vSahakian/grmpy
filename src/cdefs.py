######Class definition######
#VJS 6/2016

'''
Generic classes that are useful for this package...
'''


class db:
    '''
    This class describes a set of events
    '''
    
    def __init__(self,event,sta,N,ml,mw,DA,DV,dist,vs30):
        '''
        Initiate the class by giving the event number (event), 
        station name (sta), station number (N), local mag (ml), moment mag (mw), 
        PGA (DA), PGV (DV), and source to site distance (dist)
        '''
        
        #Give these values to the db:
        self.evnum=event
        self.sta=sta
        self.stnum=N
        self.ml=ml
        self.mw=mw
        self.pga=DA
        self.pgv=DV
        self.r=dist
        self.vs30=vs30
        
    def plot_apga(self):
        '''
        Plots all PGA, regardless of M/r
        '''
        
        from matplotlib import pyplot as plt
        import numpy as np
        
        #Open figure
        plt.figure()
        
        #Plot...
        plt.scatter(self.mw,np.log10(self.pga),marker='o')
        plt.xlabel(r"$\mathbf{M}$")
        plt.ylabel(r"$\log_{10} PGA$")
        plt.title(r"PGA vs. $\mathbf{M}$ for all distances")
        
        plt.show()
        
    def plot_rpga(self,bmin,bmax,step):
        '''
        Plots all PGA, regardless of M/r
        Input:
            bmin:       Min value for bins
            bmax:       Max balue for bins
            step:       Step interval for bins
        '''
        
        from matplotlib import pyplot as plt
        import numpy as np
        
        #Get bins:
        bins=np.arange(bmin,bmax,step)
        #Get the bin index for each recording:
        dinds=np.digitize(np.floor(self.r),bins)
        
        #Sort data to plot:
        
        #Open figure
        f=plt.figure()
        
        #Define color scale -
        #colormap needs floats that go from 0 to 1, so bins must be normalized.
        #bins is not a float though, and in order to divide by the scalar it 
        #must first be converted to a float. 
        colors=plt.cm.rainbow(bins.astype(float)/bins.max())
        #FIGURE OUT THE COLORBAR PROBLEM!!
        #f.colorbar(colors)
        
        #Plot a different color for each distance bin:
        for i in bins:
            #Find which data points are in this bin:
            bind=np.where(dinds==i+1)[0]
            #Make an array of size len(bind),4 for the colors, so that these can
            #be plotted in scatter as x,y,z (mw, pga, color).  In this bin, all
            #the colors should be the same, so tile the color for this bin i and 
            #multiply it by an array of ones. (maybe I don't even need to do this?)
            clrs=np.ones((len(bind),4))*np.tile(colors[i,:],(len(bind),1))
            #plot
            f=plt.scatter(self.mw[bind],np.log10(self.pga[bind]),edgecolors=clrs,facecolors='none',lw=0.5)

        #Label the plot - Mbold on the x, log10PGA on the y, 
        plt.xlabel(r"$\mathbf{M}$")
        plt.ylabel(r"$\log_{10} PGA$")
        plt.title(r"PGA vs. $\mathbf{M}$, binned by distance")
        
        plt.show()
        
#class hdr:
#    '''
#    This class describes the header info for a db (like AB's for now)
#    '''
#    
#    def __init__(self,dmin,dmax,mag,dist,az,year,day,hour,minu,sec,msec,nevid,idep):
#        '''
#        Has header information for events...
#        '''
#        
#        #Give header values:
#        self.dmin=dmin
#        self.dmax=dmax
#        self.mag=mag
#        self.r=dist
#        self.


#class event:
#    '''
#    This class describes one event
#    '''
#    
#    def __init__(self,event,sta,N,ml,mw,DA,DV,