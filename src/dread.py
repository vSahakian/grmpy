######Data Module######
#VJS 6/2016

#Module to read and digest data of different forms, and prepare it for grmpepy

def mread(flatfile):
    '''
    Read data from Annemarie's flatfile
    VJS 6/2016
    
    Input 
        flatfile:   String with path to the anza flatfile from AB
    Output
    
    '''
    
    import scipy.io as sio
    
    datr=sio.loadmat(flatfile)
    
    devent=datr['event']
    dsta=datr['sta']
    dhdrs=datr['hdrs']
    dN
    dMl
    dMw
    dDA
    dDV
    dPGA
    dPGV
    dlogPGA
    dlogPGV
    dnewlogPGA
    dnewlogPGV
    dVs30
    
    
    
    
    