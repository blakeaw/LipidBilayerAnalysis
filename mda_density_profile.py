import numpy as np

#dictionary of elements and their valence electron counts - used for electron profile density
valence_dict = {'H':1,'C':4,'N':5,'O':6,'P':5}

def ElectronDensityProfile(trajectory,mda_selection, axis='z',nbins=100):
	dir_ind = 2
	if	axis is 'x':
		dir_ind=0
	elif axis is 'y':
		dir_ind=1
	indices = mda_selection.indices
	natoms = len(mda_selection)
	#build the charge array
	charges = np.zeros(natoms)
	a=0
	for atom in mda_selection:
		element = atom.name[0]
		valence = valence_dict[element]
		partial = atom.charge
		total = valence+partial
		charges[a]=total
		a+=1
	zpos_list = []
	charge_list = []
	nframes = len(trajectory)
	#get the maximum box dimension along axis
	bzm=0.0
	for frame in trajectory:
		bzc = frame.dimensions[dir_ind]
		if bzc>bzm:
			bzm=bzc
	#build the profile axis
	edges = np.linspace(0.0,bzm,(nbins+1),endpoint=True)
	incr = edges[1]-edges[0]
	incr_h = incr/2.0
	centers = np.zeros(nbins)
	nedges = len(edges)
	for i in xrange(1,nedges):
		j=i-1
		centers[j]=edges[j]+inc_h
		
	counts = np.zeros(nbins)
	for frame in trajectory:
		
		sel_pos = frame._pos[indices]
		zpos = sel_pos[:,dir_ind]
		zpos_list.append(zpos)
		charge_list.append(charges)

	zpos = np.vstack(zpos_list)
	chrg = np.vstack(charge_list)
	zpos = zpos.flatten()
	chrg = chrg.flatten()
	hist,bin_edges = np.histogram(zpos,bins=nbins,weights=chrg)
	
