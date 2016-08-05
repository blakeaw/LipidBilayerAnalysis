import numpy as np

#dictionary of elements and their valence electron counts - used for electron profile density
valence_dict = {'H':1,'C':4,'N':5,'O':6,'P':5}

def ElectronDensityProfile(trajectory,mda_selection, fstart=0,fend=-1,fstep=1, axis='z',nbins=100,reference=0.0,refsel=None):
	lat_ind = [0,1]
	dir_ind = 2
	if	axis is 'x':
		dir_ind=0
		lat_ind=[1,2]
	elif axis is 'y':
		dir_ind=1
		lat_ind=[0,2]
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
	if fstart<0:
		fstart+=nframes
	if fend < 0:
		fend+=nframes+1
	#get the maximum box dimension along axis
	bzm=0.0
	nframes = 0
	sel_z_avg=0.0
	for frame in trajectory[fstart:fend:fstep]:
		bzc = frame.dimensions[dir_ind]
		if bzc>bzm:
			bzm=bzc
		ref_sel_z = 0.0
		if refsel is not None:
			ref_com = refsel.center_of_mass()
			ref_sel_z = ref_com[dir_ind]
			sel_z_avg+=ref_sel_z
		nframes+=1
	if refsel is not None:
		reference=sel_z_avg/nframes
	#build the profile axis
	edges = np.linspace(0.0,(bzm),(nbins+1),endpoint=True)
	incr = edges[1]-edges[0]
	incr_h = incr/2.0
	centers = np.zeros(nbins)
	nedges = len(edges)
	for i in xrange(1,nedges):
		j=i-1
		centers[j]=edges[j]+incr_h
		
	counts = np.zeros(nbins)
	
	for frame in trajectory[fstart:fend:fstep]:
		
		bx = frame.dimensions[lat_ind[0]]
		by = frame.dimensions[lat_ind[1]]
		binvolume = incr*bx*by
		counts_f = np.zeros(nbins)
		sel_pos = frame._pos[indices]
		zpos = sel_pos[:,dir_ind]
		push_index = zpos/incr
		j=0
		for i in push_index:
			ii = int(np.floor(i))
			if ii >= nbins:
				ii = nbins-1
			elif ii < 0:
				ii = 0
			counts_f[ii]+=charges[j]
			j+=1
		counts_f/=binvolume
		counts+=counts_f
	counts/=nframes
	centers-=reference
	return counts,centers
	
