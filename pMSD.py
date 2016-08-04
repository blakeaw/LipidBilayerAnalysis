#we are going to use the MDAnalysis to read in topo and traj
#numpy
import numpy as np
#import my running stats class
from RunningStats import *
# import the coordinae wrapping function
from pUnwrap import mda_wrap_coordinates
#function to compute the mean square displacement (MSD)
# of a single MDAnalysis atom selection (atom_sel)
# returns 1d numpy array with 2 elements: [0]=msd [1]=std-dev
def MSD_single (trajectory, atom_sel):
	#get the number of frames from the trajectory
	nframes = len(trajectory)
	#initialize a numpy array to hold the centers of mass at each frame
	comlist = np.zeros((nframes, 3))
	#index counter
	comit = 0
	#loop over the frames and store the com |r| value
	for ts in trajectory:
		#compute the center of mass for the current frame
		com = atom_sel.center_of_mass()
		#store it in the frame list
		comlist[comit]=com
		comit+=1
	#get the com for the first frame - the reference
	com0 = comlist[0]
	# compute the vector difference r(t)-r0
	dr = comlist - com0
	# square the elements of the vector difference
	dr = dr**2
	#initialize a running stats object
	drs_stat = RunningStats()
	# loop over the frames starting at index 1
	for	i in xrange(1, len(dr)):
		#get the current vector difference 
		drc=dr[i]
		#initialize a running sum		
		drs_mag = 0.0
		#loop over the elements of the vector difference and sum them--i.e. to get |r(t) - r0|**2
		for y in drc:
			drs_mag+=y
		# push the current |r(t) - r0|**2 to the running stats object (for averaging)
		drs_stat.Push(drs_mag)
	# get the msd value 
	msd = np.zeros(2)
	msd[0] = drs_stat.Mean()
	msd[1] = drs_stat.Deviation()
	#return the msd
	return msd

#function to compute the mean square displacement (MSD)
# of a list of MDAnalysis atom selections (atom_sel_list)
# returns 2d numpy array with len(atom_sel_list)X2 elements: [*,0]=msd [*,1]=std-dev
def MSD_list (trajectory, atom_sel_list):
	#get the number of frames from the trajectory
	nframes = len(trajectory)
	#get the number of atomselections
	nsels = len(atom_sel_list)
	#initialize a numpy array to hold the center of mass vectors
	comlist = np.zeros((nsels, nframes, 3))
	#index counter for the frame number
	comit = 0
	# loop over the trajectory
	for ts in trajectory:
		#loop over the selections
		for i in xrange(nsels):
			#compute the center of mass of the current selection and current frame
			com = atom_sel_list[i].center_of_mass()
			#add to the numpy array
			comlist[i,comit]=com
		comit+=1
	#initialize a numpy array to hold the msd for each selection		
	msd = np.zeros((nsels, 2))
	#initialize a running stats object to do the averaging
	drs_stat = RunningStats()
	#loop over the selections
	for	i in xrange(nsels):
		# get the current com frame list
		comlcurr = comlist[i]
		#compute the dr list
		dr = comlcurr-comlcurr[0]
		# dr**2
		dr*=dr		
		#loop over frames starting at index 1
		for	j in xrange(1, len(dr)):
			#running sum
			drs_mag = 0.0
			#get the current dr**2 vector
			drc = dr[j]
			#loop over dr**2 vector elements and sum--i.e. to get |r(t) - r0|**2
			for y in drc:
				drs_mag+=y
			#push to the running stat object
			drs_stat.Push(drs_mag)
		#get the msd for the current selection
		msdcurr = drs_stat.Mean()
		devcurr = drs_stat.Deviation()
		#push to the msd array
		msd[i,0]=msdcurr
		msd[i,1]=devcurr
		#reset the running stats object--prepare for next selection
		drs_stat.Reset()
	#return msd array
	return msd 

#function to compute the mean square displacement (MSD)
# of a list of MDAnalysis atom selections (atom_sel_list)
# returns 2d numpy array with len(atom_sel_list)X2 elements: [*,0]=msd [*,1]=std-dev
def MSD_list_unwrap_all (trajectory, atom_sel_list, unwrap=True, verbose=False):
	#get the number of frames from the trajectory
	nframes = len(trajectory)
	#get the number of atomselections
	nsels = len(atom_sel_list)
	#initialize a numpy array to hold the center of mass vectors
	comlist = np.zeros((nsels, nframes, 3))
	#index counter for the frame number
	comit = 0
	natoms = trajectory.n_atoms
	oldcoord = np.zeros((natoms,3))
	firstframe = True
	# loop over the trajectory
	for ts in trajectory:
		if	verbose:
			print " "
			print "frame ",ts.frame 
	#unwrap coordinates -- currently unwraps all the coordinates
		if	unwrap:	
			if	verbose:
				print "unwrapping frame ",ts.frame 		
			currcoord = ts._pos[:]
			if firstframe:
				oldcoord = currcoord
				firstframe = False
			else:
				abc = ts.dimensions[0:3]
				wrapcoord = mda_wrap_coordinates(abc, currcoord, oldcoord)
				ts._pos[:] = wrapcoord[:]

		#loop over the selections
		for i in xrange(nsels):
			if	verbose:
				print "frame ",ts.frame," getting com of selection ",atom_sel_list[i] 
			#compute the center of mass of the current selection and current frame
			com = atom_sel_list[i].center_of_mass()
			#add to the numpy array
			comlist[i,comit]=com
		comit+=1
	#initialize a numpy array to hold the msd for each selection		
	msd = np.zeros((nsels, 2))
	#initialize a running stats object to do the averaging
	drs_stat = RunningStats()
	#loop over the selections
	for	i in xrange(nsels):
		# get the current com frame list
		comlcurr = comlist[i]
		#compute the dr list
		dr = comlcurr-comlcurr[0]
		# dr**2
		dr*=dr		
		#loop over frames starting at index 1
		for	j in xrange(1, len(dr)):
			#running sum
			drs_mag = 0.0
			#get the current dr**2 vector
			drc = dr[j]
			#loop over dr**2 vector elements and sum--i.e. to get |r(t) - r0|**2
			for y in drc:
				drs_mag+=y
			#push to the running stat object
			drs_stat.Push(drs_mag)
		#get the msd for the current selection
		msdcurr = drs_stat.Mean()
		devcurr = drs_stat.Deviation()
		#push to the msd array
		msd[i,0]=msdcurr
		msd[i,1]=devcurr
		if	verbose:
				print "selection number ",i," has msd ",msdcurr," with deviation ",devcurr 
		#reset the running stats object--prepare for next selection
		drs_stat.Reset()
	#return msd array
	return msd 

#function to compute the mean square displacement (MSD)
# of a list of MDAnalysis atom selections (atom_sel_list)
# returns 2d numpy array with len(atom_sel_list)X2 elements: [*,0]=msd [*,1]=std-dev
def MSD_list_unwrap_sel (trajectory, atom_sel_list, unwrap=True, verbose=False):
	#get the number of frames from the trajectory
	nframes = len(trajectory)
	#get the number of atomselections
	nsels = len(atom_sel_list)
	#initialize a numpy array to hold the center of mass vectors
	comlist = np.zeros((nsels, nframes, 3))
	#index counter for the frame number
	comit = 0

	#combine all the selections into one (for wrapping)
	msel = atom_sel_list[0]
	for s in xrange(1, nsels):
		 msel+=atom_sel_list[s]
	
	natoms = len(msel)
	oldcoord = np.zeros((natoms,3))
	index = msel.indices
	firstframe = True
	# loop over the trajectory
	for ts in trajectory:
		if	verbose:
			print " "
			print "frame ",ts.frame 
	#unwrap coordinates -- currently unwraps all the coordinates
		if	unwrap:	
			if	verbose:
				print "unwrapping frame ",ts.frame 		
			currcoord = ts._pos[index]
			if firstframe:
				oldcoord = currcoord
				firstframe = False
			else:
				abc = ts.dimensions[0:3]
				wrapcoord = mda_wrap_coordinates(abc, currcoord, oldcoord)
				ts._pos[index] = wrapcoord[:]

		#loop over the selections
		for i in xrange(nsels):
			if	verbose:
				print "frame ",ts.frame," getting com of selection ",atom_sel_list[i] 
			#compute the center of mass of the current selection and current frame
			com = atom_sel_list[i].center_of_mass()
			#add to the numpy array
			comlist[i,comit]=com
		comit+=1
	#initialize a numpy array to hold the msd for each selection		
	msd = np.zeros((nsels, 2))
	#initialize a running stats object to do the averaging
	drs_stat = RunningStats()
	#loop over the selections
	for	i in xrange(nsels):
		# get the current com frame list
		comlcurr = comlist[i]
		#compute the dr list
		dr = comlcurr-comlcurr[0]
		# dr**2
		dr*=dr		
		#loop over frames starting at index 1
		for	j in xrange(1, len(dr)):
			#running sum
			drs_mag = 0.0
			#get the current dr**2 vector
			drc = dr[j]
			#loop over dr**2 vector elements and sum--i.e. to get |r(t) - r0|**2
			for y in drc:
				drs_mag+=y
			#push to the running stat object
			drs_stat.Push(drs_mag)
		#get the msd for the current selection
		msdcurr = drs_stat.Mean()
		devcurr = drs_stat.Deviation()
		#push to the msd array
		msd[i,0]=msdcurr
		msd[i,1]=devcurr
		if	verbose:
				print "selection number ",i," has msd ",msdcurr," with deviation ",devcurr 
		#reset the running stats object--prepare for next selection
		drs_stat.Reset()
	#return msd array
	return msd 

#function to compute the mean square displacement (MSD)
# of a list of MDAnalysis atom selections (atom_sel_list)
# returns 2d numpy array with len(atom_sel_list)X2 elements: [*,0]=msd [*,1]=std-dev
def MSD_lateral_list_unwrap_sel (trajectory, atom_sel_list, plane="xy", unwrap=True, verbose=False):

	ii=0
	jj=1
	#process the plane
	if	plane=="xy" or plane=="yx":
		ii=0
		jj=1
	if	plane=="yz" or plane=="zy":
		ii=1
		jj=2
	if	plane=="xz" or plane=="zx":
		ii=0
		jj=2
	plane_index = [ii, jj]
	#get the number of frames from the trajectory
	nframes = len(trajectory)
	#get the number of atomselections
	nsels = len(atom_sel_list)
	#initialize a numpy array to hold the center of mass vectors
	comlist = np.zeros((nsels, nframes, 3))
	#index counter for the frame number
	comit = 0

	#combine all the selections into one (for wrapping)
	msel = atom_sel_list[0]
	for s in xrange(1, nsels):
		 msel+=atom_sel_list[s]
	
	natoms = len(msel)
	oldcoord = np.zeros((natoms,2))
	index = msel.indices
	firstframe = True
	# loop over the trajectory
	for ts in trajectory:
		if	verbose:
			print " "
			print "frame ",ts.frame 
	#unwrap coordinates -- currently unwraps all the coordinates
		if	unwrap:	
			if	verbose:
				print "unwrapping frame ",ts.frame 		
			currcoord = ts._pos[index]
			if firstframe:
				oldcoord = currcoord
				firstframe = False
			else:
				abc = ts.dimensions[0:3]
				wrapcoord = mda_wrap_coordinates(abc, currcoord, oldcoord)
				ts._pos[index] = wrapcoord[:]

		#loop over the selections
		for i in xrange(nsels):
			if	verbose:
				print "frame ",ts.frame," getting com of selection ",atom_sel_list[i] 
			#compute the center of mass of the current selection and current frame
			com = atom_sel_list[i].center_of_mass()
			#add to the numpy array
			comlist[i,comit]=com
		comit+=1
	#initialize a numpy array to hold the msd for each selection		
	msd = np.zeros((nsels, 2))
	#initialize a running stats object to do the averaging
	drs_stat = RunningStats()
	#loop over the selections
	for	i in xrange(nsels):
		# get the current com frame list
		comlcurr = comlist[i]
		#compute the dr list
		dr = comlcurr[:,plane_index]-comlcurr[0,plane_index]
		#print dr
		# dr**2
		dr*=dr	
		print "dr*dr", dr	
		#loop over frames starting at index 1
		for	j in xrange(1, len(dr)):
			#running sum
			drs_mag = 0.0
			#get the current dr**2 vector
			drc = dr[j,:]
			print "drc ",drc
			#loop over dr**2 vector elements and sum--i.e. to get |r(t) - r0|**2
			for y in drc:
				print "y ",y
				drs_mag+=y
			#push to the running stat object
			print "drs_mag ",drs_mag
			drs_stat.Push(drs_mag)
		#get the msd for the current selection
		msdcurr = drs_stat.Mean()
		devcurr = drs_stat.Deviation()
		#push to the msd array
		msd[i,0]=msdcurr
		msd[i,1]=devcurr
		if	verbose:
				print "selection number ",i," has msd ",msdcurr," with deviation ",devcurr 
		#reset the running stats object--prepare for next selection
		drs_stat.Reset()
	#return msd array
	return msd 
