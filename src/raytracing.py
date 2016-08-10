######Raytracing Module######
#VJS 8/2016

##Module to deal with raytracing stuff...
#Make the input files, plot results, etc.

def write_sourcein(home,run_name,veltype):
    '''
    Write a source.in file for fm3d raytracing
    VJS 8/2016
    Input:
        home:           String to home path for the main working directory
        run_name:       String to the name of the run
        veltype:        0=Vp, 1=Vs
    Output:
        source_veltype.in:     Writes out the source.in file for raytracing
    '''
    
    from os import path
    import cPickle as pickle
    import dread
    from numpy import str
    
    #Get the velocity type:
    if veltype==0:
        vm='Vp'
    elif veltype==1:
        vm='Vs'
    
    #Source is always local, for now:
    telesource=0
    
    #Get the path to the run directory:
    run_dir=path.expanduser(home+run_name+'/')
    #Get the event object directory:
    eodir=run_dir+'event_objs/'
    #Get the path to the list of event objects:
    event_list_path=eodir+run_name+'.pckl'
    
    #Also get the path for the output sources.in file:
    raydir=run_dir+'raytracing/'
    sourcefile=raydir+'sources_'+vm+'.in'
    
    ##
    #Read in the list of event objects:
    eobjs=dread.read_obj_list(event_list_path)

    f=open(sourcefile,'w')
    
    #Get the number of sources:
    numsources=len(eobjs)
    #Get line to write:
    sline=str(numsources)+'\t\t\t\t number of sources\n'
    #Write to file:
    f.write(sline)
    
    #Loop over the event objects, and write each source (each event):
    for event_ind in range(len(eobjs)):
        #Get the event:
        eventi=eobjs[event_ind]
        
        #Get the source info:
        #position
        sdepth=eventi.edepth[0]
        slat=eventi.elat[0]
        slon=eventi.elon[0]
        #number of paths
        npaths=len(eventi.sta)
        
        #Number of sections on the path:  Always 1 here, since one layer
        nsections=1
        #Define path sections:
        psections=[0,2]
        #Velocity type:  #Specified as input (veltype)
        
        #Get lines to write:
        teleline=str(telesource)+'\t\t\t\t source is local/teleseismic (0/1)\n'
        posline=str(sdepth)+' '+str(slat)+' '+str(slon)+'\t\t\t position depth(km),lat(deg),lon(deg)\n'
        npathsline=str(npaths)+'\t\t\t\t number of paths from this source\n'
        nsectionsline=str(nsections)+'\t\t\t\t number of sections on the path\n'
        psectionsline=str(psections[0])+' '+str(psections[1])+'\t\t\t\t\t define the path sections\n'
        vtypeline=str(veltype)+'\t\t\t\t define the velocity type along the path\n'
        
        #Write them:
        f.write(teleline)
        f.write(posline)
        f.write(npathsline)
        f.write(nsectionsline)
        f.write(psectionsline)
        f.write(vtypeline)
        
    f.close()
        
    #Print statement
    print 'Wrote file '+sourcefile+' for velocity type '+vm


def write_receiverin():
    '''
    Write a receiver.in in file for fm3d raytracing
    Input:
        home:           String to the home path for the working directory
        run_name:       String for the run name (db/model combo)
    Output:
        receivers.in:   Receivers.in file for fm3d raytracing
    '''
    
    from os import path
    import cPickle as pickle
    import dread
    from numpy import str
    
    #Get the path to the run directory:
    run_dir=path.expanduser(home+run_name+'/')
    #Get the event object directory:
    sodir=run_dir+'event_objs/'
    #Get the path to the list of event objects:
    station_list_path=sodir+run_name+'.pckl'
    
    #Also get the path for the output sources.in file:
    raydir=run_dir+'raytracing/'
    sourcefile=raydir+'receivers.in'
    
    
    
###################################
#Function to read in a rayfile:

def parse_rayfile(rayfile):
    '''
    Parse rays.dat output file from fm3d
    Input:
        rayfile:        Path to the rays.dat file
    Output:
        path_list:      List of arrays, each array has path information
        receiver_id:    Array with id's of the receiver from receiver.in file
        source_id:      Array with id's of the source for that path from source.in file
    '''
    
    from numpy import array,r_,expand_dims
    
    #Load file
    f=open(rayfile,'r')
    
    #intialize things
    receiver_id=array([])
    source_id=array([])
    path_list=[]
    
    #Loop until end
    while True:
        
        #Read ray headers, this line has source and station IDs
        line=f.readline()
        
        #Have I recahed the end of the file?
        if line=='':
            break
            
        receiver_id=r_[receiver_id,int(line.split()[0])]
        source_id=r_[source_id,int(line.split()[1])]
        
        #Next line has number of elements in the path
        line=f.readline()
        Nelements=int(line.split()[0])
        
        #Initalize path coords variable
        path_coordinates=array([]) 
        
        #read path elements
        for kelement in range(Nelements):
            line=f.readline()
            element=array([float(line.split()[0]),float(line.split()[1]),float(line.split()[2])])
            #turn to row vector for appending
            element=expand_dims(element,0)
            if kelement==0:
                path_coordinates=element
            else:
                path_coordinates=r_[path_coordinates,element]
            
        #Append to stupid list
        path_list.append(path_coordinates)

    f.close()

    return path_list,receiver_id,source_id
        
        

            