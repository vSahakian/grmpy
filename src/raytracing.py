######Raytracing Module######
#VJS 8/2016

##Module to deal with raytracing stuff...
#Make the input files, plot results, etc.

def write_receiverin():
    '''
    Write a test receiver in file
    '''
    
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
        
        
        