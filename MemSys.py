"""
This MemSys module defines various classes and functions used to process and analyze a lipid 
bilayer trajectory. This module assumes the structure and trajectory are initiallaly stored in MDAnalysis
objects and therefore processes MDAnalysis objects. The lipids constituting the bilayer are read in from
the MDAnalysis objects and are converted to center of mass (COM) representations. Lipids are partitioned
into an 'upper' and a 'lower' leaflet based on the z-position of the COM. The built-in analysis functions
then operate on the COM representations to compute quantities such as the lateral mean squared displacement.  
Many analysis functions allow specification of the leaflet and type of lipid to perform the its analysis on.
The primary (parent class) is the MemSys class. The analysis functions are members of the MemSys class.

Example initialization:
    import MemSys as ms
    mem_sys = ms.MemSys(mda_universe.trajectory,mda_selection_of_bilayer_lipids)
 
"""

#imports
import numpy as np
import matplotlib.cm as cm
import os
import sys
import shutil
import shelve
import multiprocessing as mp
from scipy.spatial import Voronoi
from scipy.spatial import Delaunay
#import copy

#import my running stats class
from RunningStats import *
# import the coordinate wrapping function--for unwrapping
from pUnwrap import mda_wrap_coordinates,mda_wrap_coordinates_parallel

# assumes that a 1d numpy array of floats is pass as input, but 
# does not check this
def GenRunningAverage(onednparray):
    """    
        Genates a running average array corresponding to  
        the data in a 1d numpy array.

        Parameters
        ----------
        onednparray : a 1d numpy array, assumed to be array of floats     
        
        Returns
        -------
        2d numpy array of dim len(onednparray)x2
            2dnparray[i][0] = running mean at i
            2dnparray[i][1] = running standard deviation at i    
                {i = 0; i < len(onednparray)}
    """
    averager = RunningStats()
    nele = len(onednparray)
    output = np.zeros((nele,2))
    for i in xrange(nele):
        averager.Push(onednparray[i])
        run_avg = averager.Mean()
        run_dev = averager.Deviation()
        output[i,0] = run_avg
        output[i,1] = run_dev
    return output

# This function is incomplete!
def ColorizeStepVectorClusters(vectors):
    nvecs = len(vectors)
    np.zeros(nvecs,dtype=np.int)
    colors_out = np.zeros(nvecs)
    return "nothing yet!"
    
    
class LipidCOM:
    """    
        A lipid center of mass (COM) object. This object stores the COM coordinats
        of a lipid (or other molecule or group of atoms) computed from both the wrapped
        and unwrapped atomic coordinates. This object also stores information about the
        type of lipid as well as the total mass of the lipid.  
                
    """
    def __init__(self):
         """    
            This is the initialization function of the center of the LipidCOM object. 
            This function initializes all the LipidCOM instance attributes and assigns 
            some default values.

            Parameters
            ----------
            void
                           
            Returns
            -------
            void
        """
        # lipid type/resname or other name    
        self.type="UNK"
        # wrapped coordinates
        self.com=np.zeros(3)
        # unwrapped coordinates
        self.com_unwrap=np.zeros(3)
        # total mass
        self.mass=1.0
        return
    # The name of this function could be changed to be more desriptive, e.g. 
    # extract_com_mda_residue  
    def extract(self, mda_residue, unwrap=False):
        """    
            This function "extracts" the center of mass (COM) of an MDAnalysis residue. 
            This function calls the MDAnalysis member function center_of_mass() of the residue
            to compute the center of mass of the atoms constituting the residue.

            Parameters
            ----------
            mda_residue : an MDAnalysis residue object
            unwrap : bool, Optional
                      False (default) - The COM coordinates are stored in the 
                      container designated for the unwrapped coordinate representation.
                      True - The COM coordinates are stored in the container designated
                      for the wrapped coordinate representation
                           
            Returns
            -------
            void
        """
        if unwrap:
            self.com_unwrap = mda_residue.center_of_mass()
        else:
            self.com = mda_residue.center_of_mass(pbc=True)
            self.com_unwrap = self.com[:]
        
        self.type=mda_residue.resname
        return

# a frame object 
class Frame:
    """    
        A molecular dynamics Frame object. This object stores all the LipidCOM objects
        corresponding to a specific timestep, as well as other information about that
        timestep inluding the rectangular box dimensions, simulation time.  
                
    """
    # does not check that nlipids is an int
    def __init__(self, nlipids):
        """    
            This is the initialization function of Frame object. 
            This function initializes all the Frame instance attributes and assigns 
            some default values.

            Parameters
            ----------
            nlipids : int, The number of lipids (LipidCOM objects) that this frame contains
                           
            Returns
            -------
            void
        """
        # list to store the nlipids LipidCOM objects
        self.lipidcom = []
        # box dimensions
        self.box = np.zeros(3)
        # simulation time
        self.time = np.zeros(1)
        # frame number
        self.number = np.zeros(1,dtype=np.int)
        # initialize all the LipidCOM objects
        for i in xrange(nlipids):
            self.lipidcom.append(LipidCOM())

        return
    
    def SetBox(self, box_lengths):
        """    
            This member function is used to set the box dimensions of a Frame.

            Parameters
            ----------
            box_lengths : numpy array - 1d, 3 element numpy array containing the x,y,z box sizes 
                           
            Returns
            -------
            void
        """
        self.box = box_lengths
        return

    def SetTime(self, time):
        """    
            This member function is used to set the simulation time of a Frame.

            Parameters
            ----------
            time : float, simulation time 
                           
            Returns
            -------
            void
        """
        self.time = time
        return

    def __len__(self):
            return len(self.lipidcom)

#    def COG(self,unwrapped=False):
#        cog_out = np.zeros(3)
#        for lipid in self.lipidcom:    
#            if not unwrapped:
#                cog_out+=lipid.com    
#            else:
#                cog_out+=lipid.com_unwrap
#        cog_out/=len(self)
#        return com_out
    
    def COM(self, wrapped=True):
        """    
            This member function is used to compute the overall center of mass (COM) of a Frame.
            This function uses the LipidCOM object coordinates and masses to compute the COM of
            the frame. 

            Parameters
            ----------
            unwrap : bool, Optional
                      True (default) - The wrapped LipidCOM coordinates are used to compute 
                       the COM of the frame
                      False - The unwrapped LipidCOM coordinates are used to compute 
                       the COM of the frame
                           
            Returns
            -------
            frame_com : float, center of mass of the Frame
        """
        com_out = np.zeros(3)
        total_mass = 0.0
        for lipid in self.lipidcom:    
            if wrapped:
                com_out+=lipid.com*lipid.mass
                total_mass+=lipid.mass    
            else:
                com_out+=lipid.com_unwrap*lipid.mass
                total_mass+=lipid.mass    
        com_out/=total_mass
        return com_out

#frame wrapper - the name of this class may be changed. e.g. FrameShelve
class frames:
    """    
        This is a wrapper class for the Frame object that stores a set of Frame objects
        corresponding to a molecular dynamics trajectory. This class saves the Frame objects
        on disk using the shelve module and provides an interface to access instances of
        those Frames. This class defines an append function and some built-ins to allow integer indexing
        of the frames object (like an array) to add/get instances of Frame objects corresponding to that index.  
                
    """
    _type_error ="instance of object MemSys.frames only excepts instances of MemSys.Frame"
    
    def __init__(self,prefix='/tmp/',save=False):
        """    
            This is the initialization function of the frames object. 
            This function initializes all the frame instance attributes and assigns 
            some default values.

            Parameters
            ----------
            prefix : string, Optional; The location to store the "shelve"d Frame data.   
                     '/tmp/' (default) - The data is stored in the unix/linux tmp directory.
            save : bool, Optional; determine whether to delete the shelved Frame data after object deletion 
                   False (default) - the shelved Frame data is deleted upon calling __del__ 
                   True  - the shelved Frame data is not deleted when __del__ is called 
                           
            Returns
            -------
            void
        """
        self.nframes = 0
        self.pid = os.getpid()
        if prefix == 'Default':
            prefix = '/tmp/'
        
        if    prefix[-1] != '/':
            prefix = prefix +'/'
        path = prefix
        if    save:
            path = path+'mem_sys_frames'
        else:
            path = path+'.mem_sys_frames_'+str(self.pid)            
        self.path = path
        self.save = save
        if os.path.isdir(self.path):
            shutil.rmtree(self.path)
        os.mkdir(self.path, 0755)
        self.fs_name = self.path +'/shelf_frames.db' 
        self.frame_shelf = shelve.open(self.fs_name,flag="c", protocol=2)
        return

    def __del__(self):
        """    
            Non-standard implementation for the __del__ built-in.
            Closes the Frame shelve database file and deletes the shelved Frame
            data if the frames.save parameter is False

            Parameters
            ----------
            void
                           
            Returns
            -------
            void
        """
        self.frame_shelf.close()
        if not self.save:
            if os.path.isdir(self.path):
                shutil.rmtree(self.path)
        return

    def append(self,item):
        """    
            This member function allows tail append/addition like fucntionality for a Frame object. The new Frame
            is added to the shelve database with a key n_frames and the number of Frames is incremented by 1.

            Parameters
            ----------
            item : The instance of a Frame object to be appended
                           
            Returns
            -------
            void, TypeError:  Returns a TypeError if item passed for appending is not a Frame instance.
        """
        if isinstance(item, Frame):
            self.frame_shelf[str(self.nframes)] = item
            self.nframes+=1
            return
        else:
            return TypeError(self._type_error)
    

    def __getitem__(self,key):
        """    
            Non-standard implementation for the __getitem__ built-in to allow integer
            indexing of the frames object. This allows acces to the Frame objects by an 
            integer indexing key, which are stored in the shelve database files.

            Parameters
            ----------
            key : int - The index of the Frame object being called
                           
            Returns
            -------
            Frame_obj : This is an instance of the Frame object stored at index key (pulled from the shelve database)
        """
        if key < 0:
            key += self.nframes
        elif key > self.nframes:
            key = self.nframes-1
        
        return self.frame_shelf[str(key)]

    def __setitem__(self,key,item):
        """    
            Non-standard implementation for the __setitem__ built-in to allow integer
            indexing of the frames object. This allows the Frame stored at the index key to set. 

            Parameters
            ----------
            key : int - The index of where the input Frame should be stored.
            item : Frame object - This is an instance of a Frame object to be stored at index key.               
            Returns
            -------
            void, TypeError : This function returns a TypeError if the input item is not an instance of a Frame object
            
        """
        if not isinstance(item, Frame):
            return TypeError(self._type_error)
        if key < 0:
            key+=self.nframes
        elif key >= self.nframes:
            key = self.nframes
            self.nframes+=1
        self.frame_shelf[str(key)]=item
        return 

    def __len__(self):
        return self.nframes

    def __iadd__(self,item):
        """    
            Non-standard implementation for the __iadd__ built-in which allows a Frame object
            to be appended using the '+=' operator. 

            Parameters
            ----------
            item : Frame object - This is an instance of a Frame object to be appended.               
            Returns
            -------
        """
        self.append(item)
        return self
            
# the multiprocessor parallelized functions that get copies of this object
# still return:
#    Exception OSError: OSError(2, 'No such file or directory') in  ignored
# I'm not sure why, but it is marked as ignored and it doesn't seem to cause any problems with the Frame shelve
# database file.
class par_frames:
     """    
        This class is effectively used to generate read-only copies of the frames class, which can be passed 
        to functions that do parallelized computations over the number of frames.  
                
    """
    # fs_name does not actually get used, so it is deprecated and should probably be removed at some point.
    def __init__(self, nframes, fs_name, frame_shelve):
         """    
            This is the initialization function of the par_frames object. 
            This functions stors a copy of an existing Frame shelve file object. 

            Parameters
            ----------
            nframes : int - the number of Frames stored in the shelve database
            fs_name : string - the name (prefix) of the shelve database file
            frame_shelve : the shelve file object storing the Frames to be accessible by this object. 
                           
            Returns
            -------
            void
        """
            self.nframes = nframes            
            self.fs_name = fs_name 
            #print "par_frames instance"
            #print "self.nframes ",self.nframes
            #print "self.fs_name ",self.fs_name
            #self.frame_shelf = shelve.open(self.fs_name,flag="r", protocol=2)
            self.frame_shelf = frame_shelve
            return

#    def __del__(self):
#        self.frame_shelf.close()
#        return
    
    def __getitem__(self,key):
        """    
            Non-standard implementation for the __getitem__ built-in to allow integer
            indexing of the par_frames object. This allows acces to the Frame objects stored in the shelve database using 
            an integer indexing key.

            Parameters
            ----------
            key : int - The index of the Frame object being called
                           
            Returns
            -------
            Frame_obj : This is an instance of the Frame object stored at index key (pulled from the shelve database)
        """
        if key < 0:
            key += self.nframes
        elif key > self.nframes:
            key = self.nframes-1
        
        return self.frame_shelf[str(key)]

    def __len__(self):
        return self.nframes

            
# leaflet object    
class Leaflet:
    """    
        This class object is used to group lipids together according to their bilayer leaflet. It is primarily meant to 
        store the indices of LipidCOMs as they are in a Frame.lipidcom list. This class also 
        creates sub-groups within the Leaflet based on the LipidCOM.type using LipidGroup objects.
                
    """
    def __init__(self, name):
        """    
            This is the initialization function of Leaflet object. 
            This functions initializes the lists and dicts necessary to hold 
            the Leaflet data.  

            Parameters
            ----------
            name : string - the name of the bilayer leaflet being initialized ('upper' and 'lower' are used by the MemSys class)
                           
            Returns
            -------
            void
        """
        #the name of the leaflet - e.g. 'upper' or 'lower'
        self.name = name
        #initialize a list to store the indices of lipids assigned to this leaflet
        self.members = []
        #initialize a list to hold the LipidGroup objects 
        self.groups = []
        #initialize a dictionary to store the self.groups index of LipidGroup objects
        self.group_dict = {}
        return

    def __str__(self):
        return '%s leaflet of a Membrane System with %s members and %s lipid groups' % (self.name, len(self.members), len(self.groups)) 

    def __repr__(self):
        return '%s leaflet of a Membrane System with %s members and %s lipid groups' % (self.name, len(self.members), len(self.groups))
 
    def __len__(self):
        return len(self.members)

    #consider changing var name of input 'resname' to something that doesn't conflict with LipidCOM.type  
    def AddMember(self, index, resname):
        """    
            This member function allows new lipids (by Frame.lipidcom index) to be added to the Leaflet. 

            Parameters
            ----------
            index : The index of the lipid being added to the Leaflet
            resname : the resname (or LipidCOM.type) of the lipid being added.
                           
            Returns
            -------
            void
        """
        if len(self.members) == 0:
            self.members.append([index, resname])
            self.groups.append(LipidGroup(resname))
            self.groups[0].AddMember(index)
            self.group_dict.update({resname:0})
        else:
            self.members.append([index, resname])
            addgroup = True
            group_ind = 0
            for rn in self.groups:
                if resname == rn.lg_name:
                    addgroup = False
                    break
                group_ind+=1
            if addgroup:
                self.groups.append(LipidGroup(resname))
                ng = len(self.groups)
                self.groups[ng-1].AddMember(index)
                self.group_dict.update({resname: ng-1})
            else:
                self.groups[group_ind].AddMember(index)
            
            #self.members=sorted(self.members,key=lambda self.members:self.members[1])
            
        return

    def GetGroupIndices(self, group_name):
        """    
            This member function returns the list of indices grouped in the LipidGroup object
            with LipidGroup.lg_name matching the input name. This allows for selections of LipidCOMs of a specific type. 

            Parameters
            ----------
            group_name : string - The name of the group (resname of the lipids) that indices are to returned. 
                                  Passing the string 'all' will return indices of all the lipids assigned to
                                  the Leaflet instance. If the group_name is not recognised (i.e. is not in the group_dict)
                                  The function defaults to 'all'.
                           
            Returns
            -------
            void
        """
        indices = []
        if group_name == "all":
            for element in self.group_dict:
                gindex = self.group_dict[element]
                indices += self.groups[gindex].lg_members
        elif group_name in self.group_dict:
            gindex = self.group_dict[group_name]
            indices = self.groups[gindex].lg_members
        else:
            #unkwown group name- print warning and use the default "all"
            print "!! Warning - request for unknown Lipid Group \'",group_name,"\' from the ",self.name," leaflet"
            print "!! using the default \"all\""
            for element in self.group_dict:
                gindex = self.group_dict[element]
                indices += self.groups[gindex].lg_members

        return list(indices)

    def GetMemberIndices(self):
        """    
            This member function returns the list of indices for the lipids grouped in the Leaflet instance.

            Parameters
            ----------
            void
                           
            Returns
            -------
            indices : list - a list of integer indices of the lipids in the Leaflet instance
        """
        indices = []
        for element in self.members:
            indices.append(element[0])

        return list(indices)

    def HasGroup(self, group_name):
        """    
            This member function provides a way to check if there is LipidGroup in the Leaflet instance with the input
            name.

            Parameters
            ----------
            group_name : string - The name to checked against existing LipidGroup names
                           
            Returns
            -------
            answer : bool - True if there is a LipidGroup with name group_name, False otherwise
        """
        return [group_name in self.group_dict]

    def NumGroups(self):
        """    
            This member function returns the number of unique LipidGroups that have initialized within
            an instance Leaflet

            Parameters
            ----------
            none
                           
            Returns
            -------
            number_of_groups : int - The number of unique LipidGroups
        """
        return len(self.groups)

    def GetGroupNames(self):
         """    
            This member function returns the list of LipidGroup names that current exist in the 
            the Leaflet instance

            Parameters
            ----------
            void
                           
            Returns
            -------
            names : list - a list of string LipidGroup names
        """
        return [group.lg_name for group in self.groups]

        
class LipidGroup:
    """    
        This class object is used to group lipids together according to their type/resname/name. It is primarily meant to 
        store the indices of the LipidCOMs as they are in a Frame.lipidcom list. 
        Lipid members are added dynamically using the AddMember function. 
                
    """
    def __init__(self, name):
        """    
            This is the initialization function of LipidGroup object. 
            This functions initializes the list to store its members indices.
            This function also sets the name of the LipidGroup object instance. 
           

            Parameters
            ----------
            name : string - the name/type/resname of the lipids being grouped in this object 
                           
            Returns
            -------
            void
        """
        #initialize a list to hold the member indices
        self.lg_members = []
        # the name of this lipid group
        self.lg_name = name
        return

    def AddMember(self, new_mem):
         """    
            This member function allows dynamic addition (via appending to the member list) of 
            lipids via their index to the current LipidGroup instance. 

            Parameters
            ----------
            new_mem : int - the index of the lipid being added to this lipid group 
                           
            Returns
            -------
            void
        """
        self.lg_members.append(new_mem)
        return

    def name(self):
        """    
            This a member function to return the name of the current LipidGroup instance. 

            Parameters
            ----------
            void
                           
            Returns
            -------
            name : string - the name of the lipid group (i.e. lg_name) 
        """
        return self.lg_name

def MSD_frames(frames, fstart, fend, indices, refframe, plane):
     """
        This function allows the mean squared displacement (MSD) to be computed
        for a specified subset of the Frames in a frames (or par_frames) object.     
        This function was created to be called from the function MemSys.CalcMSD_parallel
        as a function to be passed to the multiprocessor threads.   
        
        Parameters
        ----------
        frames : frames or par_frames object  - object containing all the Frames of the trajectory
        fstart : int - the first frame to start the analysis on
        fend : int - the last frame to analyze
        indices : list - list of integer indices of the LipidCOMs to include in the computation
        refframe : int - the index of the frame that is to be taken as the reference for the MSD computation
        plane : list - list of the indices corresponding to the coordinate planes (x: 0,y 1,z :2) to be included in the MSD   
        
        Returns
        -------
        msd_results - numpy array (floats) - This is a num_framesx4 numpy array containing the 
                      results of the MSD computation for the specified frames
                      msd_results[i,0] = simulation time for frame f = i + fstart
                      msd_results[i,1] = the configurational average MSD over the specified LipidCOMs for frame f = i + fstart
                      msd_results[i,2] = the standard deviation of the configurational average MSD over the specified LipidCOMs for frame f = i + fstart
                      msd_results[i,3] = an estimate of the corrsponding diffusion constant based on  
                                         the configurational average MSD over the specified LipidCOMs for frame f = i + fstart
                       {i = 0; i < num_frames}
    """
    #initialize an array to hold the ouptut
    nfc = fend - fstart + 1
    output = np.zeros((nfc,4))
    # number of lipids in the selection
    n_com = len(indices)
    #initialize a running stats object to do the configuration averaging
    drs_stat = RunningStats()
    # initialize an np array to hold coordinates for the selection
    # at the reference frame
    com_ref = np.zeros((n_com,2))
    ref_frame = frames[refframe]
    count=0
    # get the coordinates
    for i in indices:
        com_i = ref_frame.lipidcom[i].com_unwrap[plane]
        com_ref[count]=com_i[:]
        count+=1
    time_ref = ref_frame.time
    #print "nframes ",len(frames)
    #print "process; fstart ",fstart," fend ",fend
    #print "process; loop range "
    #print range(fstart,(fend+1))
    # now begin loop over the frames for this process        
    for f in range(fstart, (fend+1)):
        # get the current frame
        curr_frame = frames[f]
        # get the coordinates for the selection at this frame
        com_curr = np.zeros((n_com,2))
        count=0
        for i in indices:
            com_i = curr_frame.lipidcom[i].com_unwrap[plane]
            com_curr[count]=com_i[:]
            count+=1
        #current time
        tc = curr_frame.time
        dr = com_curr - com_ref
        drs = dr*dr
        #loop over the selections for this frame
        for    val in drs:
            drs_curr = val[:]    
            drs_mag = drs_curr.sum()
            drs_stat.Push(drs_mag)
        #get the msd for the current selection
        msdcurr = drs_stat.Mean()
        devcurr = drs_stat.Deviation()
        drs_stat.Reset()
        findex = f-fstart
        output[findex,0]=tc
        output[findex,1]=msdcurr
        output[findex,2]=devcurr
        dt = tc - time_ref
        DiffCon = 0.0
        if f != 0:
            DiffCon = msdcurr/(4.0*dt)
        output[findex,3]=DiffCon
    #    print "msdcurr ",msdcurr," DiffCon ",DiffCon
    return output

#function to compute the thickness of the membrane (in the normal direction). The algorithm is based on  
# the GridMAT-MD bilayer thickness calculation (except without the gridding procedure) 
def Thickness_frames(frames, fstart, fend, leaflets, nlipids, plane, norm):
    """
        This function allows the bilayer "thickness" to be computed
        for a specified subset of the Frames in a frames (or par_frames) object.     
        This function was created to be called used in the function MemSys.CalcThickness_parallel
        as a function to be passed to the multiprocessor threads.   
        
        Parameters
        ----------
        frames : frames or par_frames object  - object containing all the Frames of the trajectory
        fstart : int - the first frame to start the analysis on
        fend : int - the last frame to analyze
        leaflets : dict - the MemSys.leaflets instance used to define the Leaflets for this calculation
                          This input should contain the two keys, 'upper' and 'lower', corresponding
                          to instances of the Leaflet class.
        nlipids : int - the total number of LipidCOMs (or lipids) in the Leaflets
        plane : list - list of the indices corresponding to the bilayer lateral coordinate planes (x: 0,y 1,z :2)
        norm : int - index corresponding to the bilayer normal coordinate plane (x: 0,y 1,z :2) 
        
        Returns
        -------
        msd_results - numpy array (floats) - This is a num_framesx4 numpy array containing the 
                      results of the MSD computation for the specified frames
                      msd_results[i,0] = simulation time for frame f = i + fstart
                      msd_results[i,1] = the configurational average MSD over the specified LipidCOMs for frame f = i + fstart
                      msd_results[i,2] = the standard deviation of the configurational average MSD over the specified LipidCOMs for frame f = i + fstart
                      msd_results[i,3] = an estimate of the corrsponding diffusion constant based on  
                                         the configurational average MSD over the specified LipidCOMs for frame f = i + fstart
                       {i = 0; i < num_frames}
    """
    #upper_match = []
    #lower_match = []
    xi = plane[0]
    yi = plane[1]
    zi = norm
    comcup = np.zeros(3)
    comclo = np.zeros(3)
    dcom = np.zeros(3)
    nfc = fend - fstart + 1
    nlc = nlipids
    zdists = np.zeros((nfc, nlc, 1))
    zmaps = np.zeros((nfc, nlc, 6))
    #dcoms = np.zeros(3)
    f=0
    times = np.zeros(nfc)
    
    for    f in range(fstart,(fend+1)):
        n=0
        fr = frames[f]
        boxc = fr.box
        boxc_xh = boxc[xi]/2.0
        boxc_yh = boxc[yi]/2.0
        dt = fr.time
        findex = f-fstart
        times[findex]=dt
        for memu in leaflets['upper'].members:
            idu = memu[0]
            comcup = fr.lipidcom[idu].com
            distxy = 10000.0
            distz = 0.0
            mindex = 0
            zlom = 0.0
            zhim = 0.0
            xavgm = 0.0
            yavgm = 0.0
            for meml in leaflets['lower'].members:
                idl = meml[0]
                comclo = fr.lipidcom[idl].com
                dcom = comcup-comclo
                dx = dcom[xi]
                dy = dcom[yi]
                dz = dcom[zi]
                #Minimum image -- coordinates must be pre-wrapped 
                if np.absolute(dx) > boxc_xh:
                    dx = boxc[xi] - np.absolute(comcup[xi]-boxc_xh) - np.absolute(comclo[xi]-boxc_xh)
                if np.absolute(dy) > boxc_yh:
                    dy = boxc[yi] - np.absolute(comcup[yi]-boxc_yh) - np.absolute(comclo[yi]-boxc_yh)
                rxy = np.sqrt(dx**2+dy**2)
                #get 4d map values
                comavg = (comcup+comclo)/2.0
                xavg = comavg[xi]
                yavg = comavg[yi]
                zlo = comclo[zi]
                zhi = comcup[zi]
                if    rxy<distxy:
                    distxy=rxy
                    distz = np.absolute(dz)
                    mindex=meml
                    xavgm = xavg
                    yavgm = yavg
                    zlom = zlo
                    zhim = zhi
                    
            #upper_match.append([mindex,distz])
            #print "n ",n," xvg ", xavgm," yvg ", yavgm
            
            zdists[findex,n]=distz
            #maps
            zmaps[findex,n,0]=dt
            zmaps[findex,n,1]=xavgm
            zmaps[findex,n,2]=yavgm
            zmaps[findex,n,3]=zlom
            zmaps[findex,n,4]=zhim
            zmaps[findex,n,5]=distz
            
            n+=1
        for meml in leaflets['lower'].members:
            idl = meml[0]
            comclo = fr.lipidcom[idl].com
            distxy = 10000.0
            distz = 0.0
            mindex = 0
            zlom = 0.0
            zhim = 0.0
            xavgm = 0.0
            yavgm = 0.0
            for memu in leaflets['upper'].members:
                idu = memu[0]
                comcup = fr.lipidcom[idu].com
                dcom = comclo-comcup
                dx = dcom[xi]
                dy = dcom[yi]
                dz = dcom[zi]
                #Minimum image -- coordinates must be pre-wrapped 
                if np.absolute(dx) > boxc_xh:
                    dx = boxc[xi] - np.absolute(comclo[xi]-boxc_xh) - np.absolute(comcup[xi]-boxc_xh)
                if np.absolute(dy) > boxc_yh:
                    dy = boxc[yi] - np.absolute(comclo[yi]-boxc_yh) - np.absolute(comcup[yi]-boxc_yh)
                rxy = np.sqrt(dx**2+dy**2)
                #get 4d map values
                comavg = (comcup+comclo)/2.0
                xavg = comavg[xi]
                yavg = comavg[yi]
                zlo = comclo[zi]
                zhi = comcup[zi]
                if    rxy<distxy:
                    distxy=rxy
                    distz = np.absolute(dz)
                    mindex=meml
                    xavgm = xavg
                    yavgm = yavg
                    zlom = zlo
                    zhim = zhi
            #upper_match.append([mindex,distz])
            #print "n ",n," xvg ", xavgm," yvg ", yavgm
            zdists[findex,n]=distz
            #maps
            zmaps[findex,n,0]=dt
            zmaps[findex,n,1]=xavgm
            zmaps[findex,n,2]=yavgm
            zmaps[findex,n,3]=zlom
            zmaps[findex,n,4]=zhim
            zmaps[findex,n,5]=distz
            n+=1
        
        #break
    zavgs = np.zeros((nfc, 3))
    zdtstat = RunningStats()    
    for fr in xrange(nfc):
        currtime = times[fr]
        dt = currtime 
        curr = zdists[fr,:]
        zavgcurr = curr.mean()            
        zdevcurr = curr.std()
#        zdtstat.Push(zavgcurr)
#        zdtcurr = zdtstat.Mean()
#        zdtdcurr = zdtstat.Deviation()
        zavgs[fr,0]=dt
        zavgs[fr,1]=zavgcurr
        zavgs[fr,2]=zdevcurr
#        zavgs[fr,3]=zdtcurr
#        zavgs[fr,4]=zdtdcurr
    out = [zavgs,zmaps]
    return out
    #return zavgs
    #return zmaps
        

## this is the main class - the Membrane System (MemSys) object
class MemSys:
    # pass the mda anaylis trajectory object and a selection with the membrane (i.e. w/o water and ions)
    # optional - specify the plane that the membrane is in - default is xy with normal in z
    def __init__(self, mda_traj, mem_sel, plane="xy",fskip=1,frame_path='Default',frame_save=False,nprocs=1):
        #defaults - xy plane with z normal
        ii=0
        jj=1
        kk=2    
        if    plane=="yz" or plane=="zy":
            ii=1
            jj=2
            kk=0
        if    plane=="xz" or plane=="zx":
            ii=0
            jj=2
            kk=1
        #parallelize loading -- currently just applies to unwrapping
        parallel=False
        if nprocs>1:
            parallel=True
        #store the indices of the plane directions
        self.plane = [ii, jj]
        # store the index of the normal direction
        self.norm = kk 
            
        #initialize leaflet objects
        self.leaflets = {'upper':Leaflet('upper'),'lower':Leaflet('lower')}
        self.com_leaflet = []
    
        
        #get the number of lipids (residues)
        self.nlipids=mem_sel.n_residues
        #initialize an empty cluster list - used to store the clusters built in the last call of 'CheckClustering'
        self.clusters = [] # after 'CheckClustering' is called, the outersize len(self.clusters) should equal self.nframes
        #initialize empty frame list
        #self.frame=[]
        self.frame = frames(prefix=frame_path,save=frame_save)
        #loop over the frames
        f=0
        for frame in mda_traj[::fskip]:
            print "doing frame ",frame.frame
            #add the frame object for this frame
            cframe = Frame(self.nlipids)
            # set the box dimensions and the time for this frame
            cframe.SetBox(frame.dimensions[0:3])
            cframe.SetTime(frame.time)
            #print "time ",frame.time
            cframe.number = f
            # loop over the residues (lipids) and get the centers of mass
            r=0            
            for res in mem_sel.residues:
                cframe.lipidcom[r].extract(res)
                cframe.lipidcom[r].mass = res.total_mass()
                r+=1
            #append the frame
            self.frame.append(cframe)
            f+=1
        #get the number of frames from the trajectory
        self.nframes = f
        #now we need to unwrap the coordinates
        natoms = len(mem_sel)
        oldcoord = np.zeros((natoms,3))
        currcoord = np.zeros((natoms,3))
        wrapcoord = np.zeros((natoms,3))
        index = mem_sel.indices

        firstframe = True
        # loop over the trajectory again to get unwrapped coordinates
        # unwrap the raw residue coordinates - then get the COMs
        f=0
        for frame in mda_traj[::fskip]:    
            #first we unwrapp
            print "unwrapping frame ",frame.frame
            currcoord = frame._pos[index]
            if firstframe:
                oldcoord = np.copy(currcoord)
                firstframe = False
            else:
                abc = frame.dimensions[0:3]
                if parallel:
                    wrapcoord = mda_wrap_coordinates_parallel(abc, currcoord, oldcoord,nprocs=nprocs)
                else:
                    wrapcoord = mda_wrap_coordinates(abc, currcoord, oldcoord)
                frame._pos[index] = wrapcoord[:]
                oldcoord = np.copy(wrapcoord)
            #now we need to adjust for the center of mass motion of the membrane -- for simplicity set all frames to (0,0,0)
            # to remove center of mass motion of the membrane
            mem_com = mem_sel.center_of_mass()
            frame._pos[index] -= mem_com
            r=0    
            cframe = self.frame[f]        
            for res in mem_sel.residues:
                cframe.lipidcom[r].extract(res, unwrap=True)
                r+=1
            self.frame[f]=cframe
            f+=1            

        # now we can assign the lipids to the leaflets 
        # NOTE: Lipids are only assigned to leaflets once based on the 
        #       first frame of the trajectory
        
        #first- compute the average position along the normal direction
        zstat = RunningStats()
        for lipcom in self.frame[0].lipidcom:
            zstat.Push(lipcom.com_unwrap[self.norm])
        zavg = zstat.Mean()
        # now loop over the lipids
        l = 0        
        for lipcom in self.frame[0].lipidcom:
            pos = ""
            # decide which leaflet
            if lipcom.com_unwrap[self.norm]>zavg:
                pos = 'upper'                
            elif lipcom.com_unwrap[self.norm]<zavg:
                pos = 'lower'
            #add to the chosen leaflet
            self.com_leaflet.append(pos)
            self.leaflets[pos].AddMember(l, lipcom.type)    
            l+=1
        #complete
        return
    
    def __str__(self):
        return 'Membrane System with %s frames and %s lipids/components' % (self.nframes, self.nlipids) 
    def __repr__(self):
        return 'Membrane System with %s frames and %s lipids/components' % (self.nframes, self.nlipids)

    def NumberOfUniqueGroups(self):
        resnames = []
        for leaflet in self.leaflets:
            for group in leaflet.groups:
                gname = group.name()
                if gname not in resnames:
                    resnames.append(gname)
        return len(resnames)
    #def LeafletCOM(leaflet_name,frame_num):
        

    # function to compute the mean squared displace (msd) along with the diffusion constant of a group 
    def CalcMSD(self, leaflet="both",group="all"):
        indices = []
        #diffusion dimension - assume lateral so, dim=2
        dim=2
        if leaflet == "both":
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.GetGroupIndices(group)
        elif leaflet == "upper":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.GetGroupIndices(group)
        elif leaflet == "lower":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.GetGroupIndices(group)
        else:
            #unknown option--use default "both"
            print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the ",self.name," leaflet"
            print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.GetGroupIndices(group)
        n_com = len(indices)
        #store the coordinates of the selected LipidCOMs in a single numpy array
        selcoords = np.zeros((self.nframes,n_com,2))
        
        for f in xrange(self.nframes): 
            count=0
            for i in indices:
                com_curr = self.frame[f].lipidcom[i].com_unwrap[self.plane]
                selcoords[f,count]=com_curr[:]
                count+=1
        
        #initialize a numpy array to hold the msd for the selection        
        msd = np.zeros((self.nframes, 7))
        #initialize a running stats object to do the averaging
        drs_stat = RunningStats()
        #initialize a running stats object for the diffusion constant (frame/time average)
        diff_stat = RunningStats()
        #running stats object for time averaging
        msd_stat = RunningStats()
        #loop over the frames starting at index 1
        #print comlist
        #print len(comlist)
        coml0 = selcoords[0,:,:]
        t0 = self.frame[0].time
        #print coml0
        for i in xrange(1, self.nframes):
            # get the current com frame list
            tc = self.frame[i].time
            dt = tc
            comlcurr = selcoords[i,:,:]
            dr = comlcurr - coml0
            drs = dr*dr
            #loop over the selections for this frame
            for    val in drs:
                drs_curr = val[:]    
                drs_mag = drs_curr.sum()
                drs_stat.Push(drs_mag)
            #get the msd for the current selection
            msdcurr = drs_stat.Mean()
            devcurr = drs_stat.Deviation()
            drs_stat.Reset()
            msd_stat.Push(msdcurr)
            msd_tavg = msd_stat.Mean()
            msd_dev = msd_stat.Deviation()            
            #dt = times[i]-times[0]
            DiffCon = msd_tavg/(2.0*dim*dt)
            diff_stat.Push(DiffCon)
            #print "msdcurr ",msdcurr
            #push to the msd array
            
            msd[i,0]=dt
            msd[i,1]=msdcurr
            msd[i,2]=msd_tavg
            msd[i,3]=msd_dev
            msd[i,4]=DiffCon
            msd[i,5]=diff_stat.Mean()
            msd[i,6]=diff_stat.Deviation()
        #return msd array
        return msd 

    #function to compute the thickness of the membrane (in the normal direction). The algorithm is based on  
    # the GridMAT-MD bilayer thickness calculation (except without the gridding procedure) 
    def CalcMembraneThickness(self):
        #upper_match = []
        #lower_match = []
        xi = self.plane[0]
        yi = self.plane[1]
        zi = self.norm
        comcup = np.zeros(3)
        comclo = np.zeros(3)
        dcom = np.zeros(3)
        zdists = np.zeros((self.nframes, self.nlipids, 1))
        zmaps = np.zeros((self.nframes, self.nlipids, 6))
        #dcoms = np.zeros(3)
        f=0
        
        for f in xrange(self.nframes): 
            n=0
            fr = self.frame[f]
            boxc = fr.box
            boxc_xh = boxc[xi]/2.0
            boxc_yh = boxc[yi]/2.0
            dt = fr.time
            for memu in self.leaflets['upper'].members:
                idu = memu[0]
                comcup = fr.lipidcom[idu].com
                distxy = 10000.0
                distz = 0.0
                mindex = 0
                zlom = 0.0
                zhim = 0.0
                xavgm = 0.0
                yavgm = 0.0
                for meml in self.leaflets['lower'].members:
                    idl = meml[0]
                    comclo = fr.lipidcom[idl].com
                    dcom = comcup-comclo
                    dx = dcom[xi]
                    dy = dcom[yi]
                    dz = dcom[zi]
                    #Minimum image -- coordinates must be pre-wrapped 
                    if np.absolute(dx) > boxc_xh:
                        dx = boxc[xi] - np.absolute(comcup[xi]-boxc_xh) - np.absolute(comclo[xi]-boxc_xh)
                    if np.absolute(dy) > boxc_yh:
                        dy = boxc[yi] - np.absolute(comcup[yi]-boxc_yh) - np.absolute(comclo[yi]-boxc_yh)
                    rxy = np.sqrt(dx**2+dy**2)
                    #get 4d map values
                    comavg = (comcup+comclo)/2.0
                    xavg = comavg[xi]
                    yavg = comavg[yi]
                    zlo = comclo[zi]
                    zhi = comcup[zi]
                    if    rxy<distxy:
                        distxy=rxy
                        distz = np.absolute(dz)
                        mindex=meml
                        xavgm = xavg
                        yavgm = yavg
                        zlom = zlo
                        zhim = zhi
                        
                #upper_match.append([mindex,distz])
                #print "n ",n," xvg ", xavgm," yvg ", yavgm
                zdists[f,n]=distz
                #maps
                zmaps[f,n,0]=dt
                zmaps[f,n,1]=xavgm
                zmaps[f,n,2]=yavgm
                zmaps[f,n,3]=zlom
                zmaps[f,n,4]=zhim
                zmaps[f,n,5]=distz
                
                n+=1
            for meml in self.leaflets['lower'].members:
                idl = meml[0]
                comclo = fr.lipidcom[idl].com
                distxy = 10000.0
                distz = 0.0
                mindex = 0
                zlom = 0.0
                zhim = 0.0
                xavgm = 0.0
                yavgm = 0.0
                for memu in self.leaflets['upper'].members:
                    idu = memu[0]
                    comcup = fr.lipidcom[idu].com
                    dcom = comclo-comcup
                    dx = dcom[xi]
                    dy = dcom[yi]
                    dz = dcom[zi]
                    #Minimum image -- coordinates must be pre-wrapped 
                    if np.absolute(dx) > boxc_xh:
                        dx = boxc[xi] - np.absolute(comclo[xi]-boxc_xh) - np.absolute(comcup[xi]-boxc_xh)
                    if np.absolute(dy) > boxc_yh:
                        dy = boxc[yi] - np.absolute(comclo[yi]-boxc_yh) - np.absolute(comcup[yi]-boxc_yh)
                    rxy = np.sqrt(dx**2+dy**2)
                    #get 4d map values
                    comavg = (comcup+comclo)/2.0
                    xavg = comavg[xi]
                    yavg = comavg[yi]
                    zlo = comclo[zi]
                    zhi = comcup[zi]
                    if    rxy<distxy:
                        distxy=rxy
                        distz = np.absolute(dz)
                        mindex=meml
                        xavgm = xavg
                        yavgm = yavg
                        zlom = zlo
                        zhim = zhi
                #upper_match.append([mindex,distz])
                #print "n ",n," xvg ", xavgm," yvg ", yavgm
                zdists[f,n]=distz
                #maps
                zmaps[f,n,0]=dt
                zmaps[f,n,1]=xavgm
                zmaps[f,n,2]=yavgm
                zmaps[f,n,3]=zlom
                zmaps[f,n,4]=zhim
                zmaps[f,n,5]=distz
                n+=1
            
            #break
        zavgs = np.zeros((self.nframes, 5))
        zdtstat = RunningStats()    
        for fr in xrange(self.nframes):
            currtime = self.frame[fr].time
            dt = currtime 
            curr = zdists[fr,:]
            zavgcurr = curr.mean()            
            zdevcurr = curr.std()
            zdtstat.Push(zavgcurr)
            zdtcurr = zdtstat.Mean()
            zdtdcurr = zdtstat.Deviation()
            zavgs[fr,0]=dt
            zavgs[fr,1]=zavgcurr
            zavgs[fr,2]=zdevcurr
            zavgs[fr,3]=zdtcurr
            zavgs[fr,4]=zdtdcurr

        return zavgs,zmaps
        #return zmaps

    # a simple cluster/chain analysis routine
    def CheckClustering(self, leaflet="both",group="all", dist=10.0):
        indices = []
        
        #diffusion dimension - assume lateral so, dim=2
        dim=2
        if leaflet == "both":
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.GetGroupIndices(group)
        elif leaflet == "upper":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.GetGroupIndices(group)
        elif leaflet == "lower":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.GetGroupIndices(group)
        else:
            #unknown option--use default "both"
            print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the ",self.name," leaflet"
            print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.GetGroupIndices(group)
        n_com = len(indices)

        #print "there are ",len(indices)," members"
        xi = self.plane[0]
        yi = self.plane[1]
        zi = self.norm
        #reset the system cluster list
        self.clusters = []
        # numpy array to store output for return
        outdata = np.zeros((self.nframes,13))
        #stats objects - time averages
        ncstat = RunningStats() #number of clusters
        asstat = RunningStats() # average cluster size
        misstat = RunningStats() # minimum cluster size
        masstat = RunningStats() # maximum cluster size                
        #loop over frames        
        
        for f in xrange(self.nframes):
            fr = self.frame[f]
            ctime = fr.time
            clusters = []
#            masterlistf = []
#            masterlistf += masterlist
            #rebuild the master list each frame
            masterlistf = list()
            for i in indices:
                masterlistf.append([i, False])
#            print "master ",masterlistf
            boxc=fr.box
            boxc_xh = boxc[xi]/2.0
            boxc_yh = boxc[yi]/2.0
            #print boxc
            clustind = 0
            neighborlist = []
            while len(masterlistf)>0:
                #print "master ",masterlistf
                start = masterlistf[0][0]
                masterlistf[0][1]=True
            #    print 
                # reset the neighborlist
                neighborlist = []
                #seed the neighborlist with the start
                neighborlist.append(start)
                #now loop over the neighborlist and build neighbors and neighbors of neigbors for this cluster
                i=0
                while i < len(neighborlist):
                    ele = neighborlist[i]
                    startn = ele
                    coms = fr.lipidcom[startn].com        
                    #get neighbors of the start
                    #mindex=0
                    for j in xrange(len(masterlistf)):
                    #for elem in masterlistf:
                        elem = masterlistf[j]
                        incluster = elem[1]
                    #    print "second incluster ",incluster
                        if not incluster:
                            ci = elem[0]
                            comc = fr.lipidcom[ci].com
                            #dcom = comc-coms
                            dx = comc[xi]-coms[xi]
                            dy = comc[yi]-coms[yi]
                            #rxy = np.sqrt(dx*dx+dy*dy)
                            #print dx," ",dy," ",rxy
                            
                            #Minimum image -- coordinates must be pre-wrapped 
                            if np.absolute(dx) > boxc_xh:
                                dx = boxc[xi] - np.absolute(comc[xi]-boxc_xh) - np.absolute(coms[xi]-boxc_xh)
                            if np.absolute(dy) > boxc_yh:
                                dy = boxc[yi] - np.absolute(comc[yi]-boxc_yh) - np.absolute(coms[yi]-boxc_yh)
                            rxy = np.sqrt(dx*dx+dy*dy)
                            #print "rxy ",rxy," dx ",dx," dy ",dy
                            if    rxy <= dist:
                                #print "passed! adding ",masterlistf[mindex][0]," to the neighborlist"
                                neighborlist.append(masterlistf[j][0])
                                masterlistf[j][1]=True
                        #mindex+=1
                    i+=1
                #filter the masterlistf
            #    print "neighlist", neighborlist
                masterlistf=list([v for v in masterlistf if v[1] == False])
                if len(neighborlist) > 1:
                    clusters.append([])
                    clusters[clustind]=list(neighborlist)
                    #print "clustind clusters[clustind]"
                    #print clustind, " ",clusters
                    clustind+=1
                        
            #print masterlistf
            #filter out single points
            #clusters = [v for v in clusters if len(v) > 1]
            nclusters = len(clusters)
            clsizestat = RunningStats()
            mini = 100000000
            maxi = -1000000
            for cluster in clusters:
                size = len(cluster)
                clsizestat.Push(size)
                if    size>maxi:
                    maxi=size
                if size < mini:
                    mini=size
            avgsize = clsizestat.Mean()
            #store instantaneous values
            outdata[f,0] = ctime
            outdata[f,1]= nclusters
            outdata[f,2] = avgsize
            outdata[f,3] = mini
            outdata[f,4] = maxi
            #push to the time averages
            ncstat.Push(nclusters)
            asstat.Push(avgsize)
            misstat.Push(mini)
            masstat.Push(maxi)
            #store current time averages
            outdata[f,5] = ncstat.Mean()
            outdata[f,6] = ncstat.Deviation()
            outdata[f,7] = asstat.Mean()
            outdata[f,8] = asstat.Deviation()
            outdata[f,9] = misstat.Mean()
            outdata[f,10] = misstat.Deviation()
            outdata[f,11] = masstat.Mean()
            outdata[f,12] = masstat.Deviation()
            # now add cluster list to the system storage
            self.clusters.append(list(clusters))
            #print clusters
            print "Frame ",f
            print "There are ",nclusters," clusters with an average size of ",avgsize
            print "the largest cluster was ",maxi," and the smallest was ",mini

        return outdata    
    #takes the cluster lists from self.clusters and gets the plane coordinates 
    # need to call the 'CheckClustering' function before calling this one
    def ExportClustersForPlotting(self):
        if len(self.clusters) == 0:
            print "Warning!! - call to \'ExportClustersForPlotting\' of a MemSys object with no cluster lists"
            print " ---------- the \'CheckClustering\' function needs to be called first!"
            return 
        xi = self.plane[0]
        yi = self.plane[1]
        #get the maximum number of clusters from any of the frames
        maxsize = 0
        for f in xrange(len(self.clusters)):
            nclust = len(self.clusters[f])
            if nclust>maxsize:
                maxsize=nclust
        #generate a color array
        colors = cm.rainbow(np.linspace(0, 1, maxsize))
        output = []
        for f in xrange(len(self.clusters)):
            frame_clusters = self.clusters[f]
            frame_data = []
            nclust = len(frame_clusters)
            #print len(frame_clusters)
            #print len(colors)
            c = 0
            xcoord = []
            #xm1 = []
            #xp1 = []
            ycoord = []
            #ym1 = []
            #yp1 =[]
            coord_color = []
            for cluster in frame_clusters:
                for index in cluster:
                    xc = self.frame[f].lipidcom[index].com[xi]
                    #xcm1 = self.frame[f].lipidcom[index].com[xi]-self.frame[f].box[xi]
                    #xcp1 = self.frame[f].lipidcom[index].com[xi]+self.frame[f].box[xi]
                    yc = self.frame[f].lipidcom[index].com[yi]
                    #ycm1 = self.frame[f].lipidcom[index].com[yi]-self.frame[f].box[yi]
                    #ycp1 = self.frame[f].lipidcom[index].com[yi]+self.frame[f].box[yi]
                     xcoord.append(xc)
                    #xm1.append(xcm1)
                    #xp1.append(xcp1)
                    ycoord.append(yc)
                    #ym1.append(ycm1)
                    #yp1.append(ycp1)
                    #print c," ",colors[c]
                    coord_color.append(colors[c])
                c+=1    
            #output.append([xm1,xcoord,xp1,ym1,ycoord,yp1,coord_color])
            output.append([xcoord,ycoord,coord_color])
        return output

    # function to compute an approximation of the area per lipid of a group using 
    # closest neighbor circles 
    def CalcAreaPerLipid_ClosestNeighborCircle(self, leaflet="both",group="all"):
        
        #diffusion dimension - assume lateral so, dim=2
        dim=2
        do_leaflet = []
        nlip = 0
        if leaflet == "both":
            do_leaflet.append('upper')
            do_leaflet.append('lower')
            nlip=self.nlipids
        
        elif leaflet == "upper" or leaflet == "lower":
            do_leaflet.append(leaflet)
            nlip = len(self.leaflets[leaflet])
        else:
            #unknown option--use default "both"
            print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the ",self.name," leaflet"
            print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
            
        xi = self.plane[0]
        yi = self.plane[1]
        zi = self.norm
        sub_fact = (2.0*np.pi/3.0 - np.sqrt(3.0)/2.0)
        #initialize a numpy array to hold the msd for the selection        
        areas = np.zeros((self.nframes, 5))
        #initialize a running stats object to do the averaging
        area_stat = RunningStats()
        n_leaflet = len(do_leaflet)
        #build the index lists
        indices_leaflet = {}
        all_mem_leaflet = {}
        for leaflets in do_leaflet:
            indices = list()
            curr_leaf = self.leaflets[leaflets]
            indices+=curr_leaf.GetGroupIndices(group)
            n_com = len(indices)
            all_mem = list(self.leaflets[leaflets].GetMemberIndices())
            all_mem_leaflet[leaflets] = list(all_mem)
            indices_leaflet[leaflets]=list(indices)
            
        
        #loop over the frames
        for f in xrange(self.nframes): 
            fr = self.frame[f]
            dt = fr.time
            boxc=fr.box
            boxc_xh = boxc[xi]/2.0
            boxc_yh = boxc[yi]/2.0
            lat_area = boxc_xh*boxc_yh*4.0
            if leaflet == 'both':
                lat_area*=2.0
            
            area_stat_config = RunningStats()
            #loop over the leaflets
            for leaflets in do_leaflet:
                indices = indices_leaflet[leaflets]
                all_mem = all_mem_leaflet[leaflets]
                #loop over the group indices in this leaflet
                for    index in indices:
                    comc = fr.lipidcom[index].com[:]
                    rdist_min = 10000.0
                    #loop over the COMs of non group 
                    #get all the leaflet members
                    
                    for a in all_mem:
                        #print "a ",a
                        if a != index:
                            comn = fr.lipidcom[a].com[:]
                            dx = comc[xi]-comn[xi]
                            dy = comc[yi]-comn[yi]
                        
                            #Minimum image -- coordinates must be pre-wrapped 
                            if np.absolute(dx) > boxc_xh:
                                dx = boxc[xi] - np.absolute(comc[xi]-boxc_xh) - np.absolute(comn[xi]-boxc_xh)
                            if np.absolute(dy) > boxc_yh:
                                dy = boxc[yi] - np.absolute(comc[yi]-boxc_yh) - np.absolute(comn[yi]-boxc_yh)
                            rxy = np.sqrt(dx*dx+dy*dy)
                            #print "rxy ",rxy," dx ",dx," dy ",dy
                            if    rxy < rdist_min:
                                rdist_min = rxy
                    #got the min dist, now compute area
                    #print "rdist_min ",rdist_min
                    area = np.pi*rdist_min*rdist_min - (rdist_min*rdist_min)*sub_fact
                    area_stat_config.Push(area)
            area_conf_avg = area_stat_config.Mean()
            area_stat.Push(area_conf_avg)
            area_time_run = area_stat.Mean()
            area_time_run_dev = area_stat.Deviation()
            #print "time ",dt
            areas[f][0]=dt
            areas[f][1]=area_conf_avg
            areas[f][2]=area_time_run
            areas[f][3]=area_time_run_dev
            areas[f][4]=lat_area/nlip
        return areas

    # function to compute the area per lipid using the lateral box sizes and numbers of lipids:
    def CalcAreaPerLipid_Box(self, leaflet="both"):
        
        #diffusion dimension - assume lateral so, dim=2
        dim=2
        do_leaflet = []
        nlip = 0
        if leaflet == "both":
            do_leaflet.append('upper')
            do_leaflet.append('lower')
            nlip = []
            for leaflets in do_leaflet:
                nlip.append(float(len(self.leaflets[leaflets])))
        
        elif leaflet == "upper":
            do_leaflet.append(leaflet)
            nlip = len(self.leaflets[leaflet])
        elif leaflet == "lower":
            do_leaflet.append(leaflet)
            nlip = len(self.leaflets[leaflet])
        else:
            #unknown option--use default "both"
            print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the ",self.name," leaflet"
            print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
            
        xi = self.plane[0]
        yi = self.plane[1]
        zi = self.norm
        
        #initialize a numpy array to hold the msd for the selection        
        areas = np.zeros((self.nframes, 4))
        #initialize a running stats object to do the averaging
        area_stat = RunningStats()
        n_leaflet = len(do_leaflet)
            
        
        #loop over the frames
        for f in xrange(self.nframes): 
            fr = self.frame[f]
            dt = fr.time
            boxc=fr.box
            boxc_xh = boxc[xi]/2.0
            boxc_yh = boxc[yi]/2.0
            lat_area = boxc_xh*boxc_yh*4.0
            area_per_lip = lat_area/nlip
            if leaflet == 'both':
                area_per_lip = (lat_area/2.0)*( (nlip[0]+nlip[1])/(nlip[0]*nlip[1]))
            
            area_stat.Push(area_per_lip)
            area_time_run = area_stat.Mean()
            area_time_run_dev = area_stat.Deviation()
            areas[f][0]=dt
            areas[f][1]=area_per_lip
            areas[f][2]=area_time_run
            areas[f][3]=area_time_run_dev
        return areas

    # do Voronoi tesselation using the COMs as generators
    def VoronoiTesselate(self, leaflet="both",group="all"):
        indices = []
        
        #diffusion dimension - assume lateral so, dim=2
        dim=2
        if leaflet == "both":
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.GetGroupIndices(group)
        elif leaflet == "upper":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.GetGroupIndices(group)
        elif leaflet == "lower":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.GetGroupIndices(group)
        else:
            #unknown option--use default "both"
            print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the ",self.name," leaflet"
            print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.GetGroupIndices(group)
        n_com = len(indices)

        #print "there are ",len(indices)," members"
        xi = self.plane[0]
        yi = self.plane[1]
        zi = self.norm
        out_tess = []
        for    f in xrange(self.nframes):
            # get the current frame
            curr_frame = self.frame[f]
            # get the coordinates for the selection at this frame
            com_curr = np.zeros((n_com,2))
            count=0
            for i in indices:
                com_i = curr_frame.lipidcom[i].com_unwrap[self.plane]
                com_curr[count]=com_i[:]
                count+=1
            vor = Voronoi(com_curr)
            #out_tess.append([com_curr[:,0],com_curr[:,1],vor])
            out_tess.append(vor)
        return out_tess

    # do Delauny tesselation using the COMs as generators
    def DelaunayTesselate(self, leaflet="both",group="all"):
        indices = []
        
        #diffusion dimension - assume lateral so, dim=2
        dim=2
        if leaflet == "both":
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.GetGroupIndices(group)
        elif leaflet == "upper":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.GetGroupIndices(group)
        elif leaflet == "lower":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.GetGroupIndices(group)
        else:
            #unknown option--use default "both"
            print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the ",self.name," leaflet"
            print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.GetGroupIndices(group)
        n_com = len(indices)

        #print "there are ",len(indices)," members"
        xi = self.plane[0]
        yi = self.plane[1]
        zi = self.norm
        out_tess = []
        for    f in xrange(self.nframes):
            # get the current frame
            curr_frame = self.frame[f]
            # get the coordinates for the selection at this frame
            com_curr = np.zeros((n_com,2))
            count=0
            for i in indices:
                com_i = curr_frame.lipidcom[i].com_unwrap[self.plane]
                com_curr[count]=com_i[:]
                count+=1
            tri = Delaunay(com_curr)
            out_tess.append([com_curr[:,0],com_curr[:,1],tri])
        return out_tess


    # generate the step vectors of the center of mass--in the lateral dimensions
    def StepVector(self, leaflet="both",group="all",fstart=0,fend=-1,fstep=1000,wrapped=False):
        indices = []
        if fstart<0:
            fstart+=self.nframes
        if fend < 0:
            fend+=self.nframes

        #diffusion dimension - assume lateral so, dim=2
        dim=2
        if leaflet == "both":
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.GetGroupIndices(group)
        elif leaflet == "upper":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.GetGroupIndices(group)
    
        elif leaflet == "lower":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.GetGroupIndices(group)
            
        else:
            #unknown option--use default "both"
            print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the ",self.name," leaflet"
            print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.GetGroupIndices(group)
        n_com = len(indices)

        #print "there are ",len(indices)," members"
        xi = self.plane[0]
        yi = self.plane[1]
        zi = self.norm
        
        vec_ends_out = []
        for    f in xrange(fstart,fend+1,fstep):
            fprev = f-fstep
            # get the current frame
            curr_frame = self.frame[f]
            prev_frame = self.frame[fprev]
            # get the coordinates for the selection at this frame
            vec_ends = np.zeros((n_com,4))
            #vec_ends = []
            count=0
            for i in indices:
                com_i = curr_frame.lipidcom[i].com_unwrap[self.plane]
                com_j = prev_frame.lipidcom[i].com_unwrap[self.plane]
                com_j_w = prev_frame.lipidcom[i].com[self.plane]
                if wrapped: 
                    vec_ends[count,0]=com_j_w[0]
                    vec_ends[count,1]=com_j_w[1]
                else:
                    vec_ends[count,0]=com_j[0]
                    vec_ends[count,1]=com_j[1]
                vec_ends[count,2]=com_i[0] - com_j[0]
                vec_ends[count,3]=com_i[1] - com_j[1]
            #    vec_ends.append([com_j[0],com_j[0],com_i[0]-com_j[0],com_i[1]-com_j[1]])
                count+=1
            vec_ends_out.append(vec_ends)
            
        return vec_ends_out

    # generate the step vectors of the center of mass
    def StepVectorColors(self, leaflet="both",group="all"):
        indices = []            
        ngroups = 1
        group_names = []
        #diffusion dimension - assume lateral so, dim=2
        dim=2
        if leaflet == "both":
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.GetGroupIndices(group)
                curr_group_names = curr_leaf.GetGroupNames()
                if group == 'all':
                    for gname in curr_group_names:
                        if gname not in group_names:
                            group_names.append(gname)
                else:
                    group_names.append(group)
        elif leaflet == "upper":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.GetGroupIndices(group)
            curr_group_names = curr_leaf.GetGroupNames()
            if group == 'all':
                for gname in curr_group_names:
                    if gname not in group_names:
                        group_names.append(gname)
            else:
                group_names.append(group)
                
        elif leaflet == "lower":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.GetGroupIndices(group)
            curr_group_names = curr_leaf.GetGroupNames()
            if group == 'all':
                for gname in curr_group_names:
                    if gname not in group_names:
                        group_names.append(gname)
            else:
                group_names.append(group)
                
        else:
            #unknown option--use default "both"
            print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the ",self.name," leaflet"
            print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.GetGroupIndices(group)
                curr_group_names = curr_leaf.GetGroupNames()
                if group == 'all':
                    for gname in curr_group_names:
                        if gname not in group_names:
                            group_names.append(gname)
                else:
                    group_names.append(group)
                
        n_com = len(indices)
        ngroups = len(group_names)
        colors = cm.rainbow(np.linspace(0, 1, ngroups))
        #build color map
        cmap = {}
        n = 0
        for name in group_names:
            cmap[name] = colors[n]
            n+=1
        #pick a frame-just use first frame    
        curr_frame = self.frame[0]
        colors_out = np.zeros( (n_com,4))
        count=0
        for i in indices:
            name_i = curr_frame.lipidcom[i].type
            colors_out[count] = cmap[name_i]
            count+=1
            
        return colors_out,cmap

    def RemoveLeafletCOMmotion(self,leaflet="both"):
        do_leaflet = []
        nlip = 0
        if leaflet == "both":
            do_leaflet.append('upper')
            do_leaflet.append('lower')
            nlip = []
            for leaflets in do_leaflet:
                nlip.append(float(len(self.leaflets[leaflets])))
        
        elif leaflet == "upper":
            do_leaflet.append(leaflet)
            nlip = len(self.leaflets[leaflet])
        elif leaflet == "lower":
            do_leaflet.append(leaflet)
            nlip = len(self.leaflets[leaflet])
        else:
            #unknown option--use default "both"
            print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the ",self.name," leaflet"
            print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
            do_leaflet.append('upper')
            do_leaflet.append('lower')
        leaf_indices = {}

        for leaf in do_leaflet:
            leaf_indices[leaf]=list(self.leaflets[leaf].GetMemberIndices())
        
        
        for f in xrange(self.nframes):
            fr = self.frame[f]
            
            for leaf in do_leaflet:
                indices=leaf_indices[leaf]
                #get the leaflet COM
                lcom = np.zeros(3)
                masst = 0.0
                for i in indices:
                    lcom+=(fr.lipidcom[i].com_unwrap*fr.lipidcom[i].mass)
                    masst+=fr.lipidcom[i].mass
                lcom/=masst
                for i in indices:
                    fr.lipidcom[i].com_unwrap-=lcom
            self.frame[f]=fr
        return

    ############### multiprocessor parallelized versions of calculation member functions

    # parallelized version of CalcMSD- using the multiprocessing module 
    def CalcMSD_parallel(self, leaflet="both",group="all",nprocs=2,timeaverage=False):            

        indices = []
        #diffusion dimension - assume lateral so, dim=2
        dim=2
        if leaflet == "both":
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.GetGroupIndices(group)
        elif leaflet == "upper":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.GetGroupIndices(group)
        elif leaflet == "lower":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.GetGroupIndices(group)
        else:
            #unknown option--use default "both"
            print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the ",self.name," leaflet"
            print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.GetGroupIndices(group)
        n_com = len(indices)

        frame_ranges = []
        total_frames = self.nframes
        frames_per_proc_base = total_frames/nprocs
        left_over = total_frames % (frames_per_proc_base * nprocs)
        print "total frames ",total_frames
        print "frames per proc ",frames_per_proc_base
        print "left over ",left_over
        #assign base ranges
        for i in xrange(nprocs):
            fs = i*frames_per_proc_base
            fe = fs + frames_per_proc_base - 1
            frame_ranges.append([fs,fe])
        print "frame_ranges (pre-adjust):"    
        print frame_ranges
        #now adjust for leftovers - divide them "equally" over the processes
        lo = left_over
        while lo > 0:
            for i in xrange(nprocs):
                frame_ranges[i][1]+=1
                for j in xrange(i+1,nprocs):
                    frame_ranges[j][0]+=1
                    frame_ranges[j][1]+=1
                lo-=1
                if lo == 0:
                    break
        
        print "nprocs ",nprocs
        print "frame_ranges (post adjust): "
        print frame_ranges
        #initialize a numpy array to hold the msd for the selection        
        msd = np.zeros((self.nframes, 4))
        #
        msd_frames = MSD_frames
        #frames_local = getattr(self, 'frame')
        #shelf_local = shelve.open(self.frame.fs_name,flag="r", protocol=2)
        frames_local = par_frames(self.frame.nframes,self.frame.fs_name,self.frame.frame_shelf)
        #frames_local = par_frames(self.frame.nframes,self.frame.fs_name)
        #frames_local = par_frames(self.frame.nframes,self.frame.fs_name,shelf_local)
        plane_local = self.plane
        #create process pool
        pool = mp.Pool(processes=nprocs)
        results = [pool.apply_async(msd_frames,args=(frames_local,frame_ranges[i][0],frame_ranges[i][1],indices,0,plane_local)) for i in range(0,nprocs)]
    #    print "results:"
    #    print results
        results_ordered = [p.get() for p in results]
    #    print "results ordered: "
    #    print results_ordered
#        #collect results  into single array for return
        i = 0
    #    print "len(results_ordered) ",len(results_ordered)
        for p in results_ordered:
            fs = frame_ranges[i][0]
            fe = frame_ranges[i][1]
            #print fs, fe
            #print msd[fs:(fe+1)].shape
            #print p[:].shape
            msd[fs:(fe+1)] = p[:]
            i+=1
        pool.close()
        pool.join()
        #initialize a numpy array to hold the msd for the selection        
        msd_tavg = msd[:]
        if timeaverage:
            #regenerate the container
            msd_tavg = np.zeros((self.nframes, 6))
            # get the running time average
            tavg_msd = GenRunningAverage(msd[:,1])
            #slice together the values
            msd_tavg[:,0:4]=msd[:,:]
            msd_tavg[:,4:6]=tavg_msd[:,:]
            
            
        #shelf_local.close()
        return msd_tavg
    #function to compute the thickness of the membrane (in the normal direction). The algorithm is based on  
    # the GridMAT-MD bilayer thickness calculation (except without the gridding procedure) 
    def CalcMembraneThickness_parallel(self,nprocs=2,timeaverage=True):
        nlip = self.nlipids
        comcup = np.zeros(3)
        comclo = np.zeros(3)
        dcom = np.zeros(3)
        zdists = np.zeros((self.nframes, 3))
        zmaps = np.zeros((self.nframes, self.nlipids, 6))
        frame_ranges = []
        total_frames = self.nframes
        frames_per_proc_base = total_frames/nprocs
        left_over = total_frames % (frames_per_proc_base * nprocs)
        print "total frames ",total_frames
        print "frames per proc ",frames_per_proc_base
        print "left over ",left_over
        #assign base ranges
        for i in xrange(nprocs):
            fs = i*frames_per_proc_base
            fe = fs + frames_per_proc_base - 1
            frame_ranges.append([fs,fe])
        print "frame_ranges (pre-adjust):"    
        print frame_ranges
        #now adjust for leftovers - divide them "equally" over the processes
        lo = left_over
        while lo > 0:
            for i in xrange(nprocs):
                frame_ranges[i][1]+=1
                for j in xrange(i+1,nprocs):
                    frame_ranges[j][0]+=1
                    frame_ranges[j][1]+=1
                lo-=1
                if lo == 0:
                    break
        
        print "nprocs ",nprocs
        print "frame_ranges (post adjust): "
        print frame_ranges
        
        thick_frames = Thickness_frames
        frames_local = par_frames(self.frame.nframes,self.frame.fs_name,self.frame.frame_shelf)
        plane_local = self.plane
        norm_local = self.norm
        #create process pool
        pool = mp.Pool(processes=nprocs)
        results = [pool.apply_async(thick_frames,args=(frames_local,frame_ranges[i][0],frame_ranges[i][1],self.leaflets,nlip,plane_local,norm_local)) for i in range(0,nprocs)]
        print "results:"
    #    print results
        print "len(results) ",len(results)
        results_ordered = [p.get() for p in results]
        print "results ordered: "
    #    print results_ordered
#        #collect results  into single array for return
        i = 0
        #print "len(results_ordered) ",len(results_ordered)
        for p in results_ordered:
            fs = frame_ranges[i][0]
            fe = frame_ranges[i][1]
            print fs, fe
            #print msd[fs:(fe+1)].shape
            #print p[:].shape
            zdistf = p[0]
            zmapf = p[1]
            #print zdistf.shape," ",zmapf.shape
            zdists[fs:(fe+1)] = zdistf[:]
            zmaps[fs:(fe+1)] = zmapf[:]
            #zdists[fs:(fe+1)] = pg[:]
            i+=1
        pool.close()
        pool.join()
        #initialize a numpy array to hold the msd for the selection        
        zdist_tavg = zdists
        if timeaverage:
            #regenerate the container
            zdist_tavg = np.zeros((self.nframes, 5))
            # get the running time average
            tavg_dz = GenRunningAverage(zdists[:,1])
            #slice together the values
            zdist_tavg[:,0:3]=zdists[:,:]
            zdist_tavg[:,3:5]=tavg_dz[:,:]
            
            
        #shelf_local.close()
        return zdsit_tavg,zmaps
        #return zdist_tavg


