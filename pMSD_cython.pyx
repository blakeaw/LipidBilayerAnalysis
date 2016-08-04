#numpy
import numpy as np
cimport numpy as np
#cython
import cython
cimport cython
#from numpy cimport ndarray

#function to compute the mean square displacement (MSD)
# of a list of MDAnalysis atom selections (atom_sel_list)
# returns 2d numpy array with len(atom_sel_list)X2 elements: [*,0]=msd [*,1]=std-dev
@cython.boundscheck(False)
def MSD_list_cy (trajectory, atom_sel_list):
	cdef int nframes, nsels, comit, i, j
	cdef double y, drs_mag, msdcurr, devcurr
	cdef np.ndarray[np.double_t, ndim=3] comlist
	cdef np.ndarray[np.double_t, ndim=1] drs_stat
	cdef np.ndarray[np.double_t, ndim=2] comlcurr
	cdef np.ndarray[np.double_t, ndim=2] dr
	cdef np.ndarray[np.double_t, ndim=1] drc
	cdef np.ndarray[np.double_t, ndim=2] msd
	#get the number of frames from the trajectory
	nframes = len(trajectory)
	#get the number of atomselections
	nsels = len(atom_sel_list)
	#initialize a numpy array to hold the center of mass vectors
	
	comlist = np.zeros((nsels, nframes, 3), dtype=np.double)

	#index counter for the frame number
	comit = 0
	# loop over the trajectory
	for ts in trajectory:
		#loop over the selections
		#cdef int i
		for i in xrange(nsels):
			#compute the center of mass of the current selection and current frame
			com = atom_sel_list[i].center_of_mass()
			#add to the numpy array
			comlist[i,comit]=com
		comit+=1
	#initialize a numpy array to hold the msd for each selection	
		
	msd = np.zeros((nsels, 2))
	#initialize a running stats object to do the averaging
	#drs_stat = RunningStats()
	
	#loop over the selections
	#cdef int i
	for	i in xrange(nsels):
		# get the current com frame list
		comlcurr = comlist[i]
		#compute the dr list
		dr = comlcurr-comlcurr[0]
		# dr**2
		dr*=dr		
		#loop over frames starting at index 1
		#cdef int j
		
		drs_stat = np.zeros((len(dr)-1))
		for	j in xrange(1, len(dr)):
			#running sum
			drs_mag = 0.0
			#get the current dr**2 vector
			drc = dr[j]
			#loop over dr**2 vector elements and sum--i.e. to get |r(t) - r0|**2
			
			for y in drc:
				drs_mag+=y
			#push to the running stat object
			drs_stat[j-1]=drs_mag
		#get the msd for the current selection
		msdcurr = drs_stat.mean()
		devcurr = drs_stat.std()
		#push to the msd array
		msd[i,0]=msdcurr
		msd[i,1]=devcurr
		#reset the running stats object--prepare for next selection
		#drs_stat.Reset()
	#return msd array
	return msd 



