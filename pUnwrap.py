#numpy
import numpy as np
#MDAnalysis
import MDAnalysis as mda


def mda_wrap_coordinates(abc, coord, refcoord):
	wrapcoord = coord.copy()
	# box vectors (assuming rectangular)
	A = abc[0]
	B = abc[1]
	C = abc[2]
	Ah = A/2.0
	Bh = B/2.0
	Ch = C/2.0
	iA = 1.0/A
	iB = 1.0/B
	iC = 1.0/C
	
	natoms = len(coord)
	#list to hold shifts
	sAs = np.zeros(natoms)
	sBs = np.zeros(natoms)
	sCs = np.zeros(natoms)

	#compute the difference in the z coordinate
	dzs = coord[:,2] - refcoord[:,2]
	#compute the required shift
	i = 0
	for dz in dzs:
		shift = 0
		if	dz > Ch:
			shift-=1
			while (dz+shift*C)>Ch:
				shift-=1
		elif dz < -Ch:
			shift+=1
			while (dz+shift*C)<-Ch:
				shift+=1
		sCs[i]=shift
		i+=1
	#apply shift
	wrapcoord[:,2] += (sCs*C)
	#compute the difference in the y coordinate
	dys = coord[:,1] - refcoord[:,1]
	#compute the required shift
	i = 0
	for dy in dys:
		shift = 0
		if	dy > Bh:
			shift-=1
			while (dy+shift*B)>Bh:
				shift-=1
		elif dy < -Bh:
			shift+=1
			while (dy+shift*B)<-Bh:
				shift+=1
		sBs[i]=shift
		i+=1
	#apply shift
	wrapcoord[:,1] += (sBs*B)
	#compute the difference in the x coordinate
	dxs = coord[:,0] - refcoord[:,0]
	#compute the required shift
	i = 0
	for dx in dxs:
		shift = 0
		if	dx > Ah:
			shift-=1
			while (dx+shift*A)>Ah:
				shift-=1
		elif dx < -Ah:
			shift+=1
			while (dz+shift*A)<-Ah:
				shift+=1
		sAs[i]=shift
		i+=1
	#apply shift
	wrapcoord[:,0] += (sAs*A)
	return wrapcoord

def mda_unwrap(universe, out_file_name):
	frames = universe.trajectory
	print "unwrapping coordinates - ouput is ",out_file_name
	nframes=len(frames)
	# Setup writer to write aligned dcd file
	writer = mda.coordinates.DCD.DCDWriter(
		out_file_name, frames.n_atoms,
		0,
		1,
		frames.dt,
		remarks='Unwrapped trajectory')
	natoms = frames.n_atoms
	oldcoord = np.zeros((natoms,3), dtype=np.double)
	firstframe = True
	for frame in frames:
		currcoord = frame._pos[:]
		if firstframe:
			oldcoord = currcoord
			firstframe = False
			writer.write(universe.atoms)
		else:
			abc = frame.dimensions[0:3]
			wrapcoord = mda_wrap_coordinates(abc, currcoord, oldcoord)
			frame._pos[:] = wrapcoord[:]
			writer.write(universe.atoms)

#def wrap_coordinates(abc, coord, refcoord):
#	wrapcoord = coord.copy()
#	# box vectors (assuming rectangular)
#	A = abc[0]
#	B = abc[1]
#	C = abc[2]
#	Ah = A/2.0
#	Bh = B/2.0
#	Ch = C/2.0
#	iA = 1.0/A
#	iB = 1.0/B
#	iC = 1.0/C
#	
#	natoms = 1
#	#list to hold shifts
#	sAs = np.zeros(natoms)
#	sBs = np.zeros(natoms)
#	sCs = np.zeros(natoms)

#	#compute the difference in the z coordinate
#	dz = coord[2] - refcoord[2]
#	#compute the required shift
#	i = 0
#	shift = 0
#	if	dz > Ch:
#		shift-=1
#		while (dz+shift*C)>Ch:
#			shift-=1
#	elif dz < -Ch:
#		shift+=1
#		while (dz+shift*C)<-Ch:
#			shift+=1
#	sCs[i]=shift
#	#apply shift
#	wrapcoord[2] += (sCs*C)
#	#compute the difference in the y coordinate
#	dy = coord[1] - refcoord[1]
#	#compute the required shift
#	i = 0

#	shift = 0
#	if	dy > Bh:
#		shift-=1
#		while (dy+shift*B)>Bh:
#			shift-=1
#	elif dy < -Bh:
#		shift+=1
#		while (dy+shift*B)<-Bh:
#			shift+=1
#	sBs[i]=shift
#	i+=1
#	#apply shift
#	wrapcoord[1] += (sBs*B)
#	#compute the difference in the x coordinate
#	dx = coord[0] - refcoord[0]
#	#compute the required shift
#	i = 0
#	shift = 0
#	if	dx > Ah:
#		shift-=1
#		while (dx+shift*A)>Ah:
#			shift-=1
#	elif dx < -Ah:
#		shift+=1
#		while (dz+shift*A)<-Ah:
#			shift+=1
#	sAs[i]=shift
#	
#	#apply shift
#	wrapcoord[0] += (sAs*A)
#	return wrapcoord			
