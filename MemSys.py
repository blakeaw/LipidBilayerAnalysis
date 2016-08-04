#we are going to use the MDAnalysis to read in topo and traj
#numpy
import numpy as np
import matplotlib.cm as cm
import os, sys, shutil
import shelve
#import my running stats class
from RunningStats import *
# import the coordinate wrapping function--for unwrapping
from pUnwrap import mda_wrap_coordinates

#lipid center of mass object - stores the center of mass of a lipid/residue - stores both wrapped and unwrapped coordinates
class LipidCOM:

	def __init__(self):	
		self.type="UNK"
		self.com=np.zeros(3)
		self.com_unwrap=np.zeros(3)
		return
	#pass the MDAnalysis residue object - 
	# optional param 'unwrap' - False is default - implies that the residue is in wrapped coordinates
	#						  - True implies that the residue is in unwrapped coordinates	
	def extract(self, mda_residue, unwrap=False):
		if	unwrap:
			self.com_unwrap = mda_residue.center_of_mass()
		else:
			self.com = mda_residue.center_of_mass(pbc=True)
			self.com_unwrap = self.com[:]
		
		self.type=mda_residue.resname
		return
# a frame object 
class Frame:
	def __init__(self, nlipids):
		self.lipidcom = []
		self.box = np.zeros(3)
		self.time = np.zeros(1)
		for i in xrange(nlipids):
			self.lipidcom.append(LipidCOM())
		return
	def SetBox(self, box_lengths):
		self.box = box_lengths
		return
	def SetTime(self, time):
		self.time = time
		return
	def __len__(self):
			return len(self.lipidcom)
#frame wrapper
class frames:
	_type_error ="instance of object MemSys.frames only excepts instances of MemSys.Frame"
	
	def __init__(self,prefix='/tmp/'):
			self.nframes = 0
			if	prefix[-1] != '/':
				prefix = prefix +'/'
			self.path = prefix+'.mem_sys_frames'
			if os.path.isdir(self.path):
				shutil.rmtree(self.path)
			os.mkdir(self.path, 0755)
			self.fs_name = self.path +'/shelf_frames.db' 
			self.frame_shelf = shelve.open(self.fs_name,flag="c", protocol=2)
			return

	def __del__(self):
		self.frame_shelf.close()
		if os.path.isdir(self.path):
			shutil.rmtree(self.path)
		return

	def append(self,item):
		if isinstance(item, Frame):
			self.frame_shelf[str(self.nframes)] = item
			self.nframes+=1
			return
		else:
			return TypeError(self._type_error)
	

	def __getitem__(self,key):
		if key < 0:
			key += self.nframes
		elif key > self.nframes:
			key = self.nframes-1
		
		return self.frame_shelf[str(key)]

	def __setitem__(self,key,item):
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
		self.append(item)
		return self
			

# leaflet object	
class Leaflet:
	def __init__(self, name):
		self.name=name
		self.members = []
		#self.group_names = []
		#self.groups = {}
		self.groups = []
		self.group_dict = {}
		return
	def __str__(self):
		return '%s leaflet of a Membrane System with %s members and %s lipid groups' % (self.name, len(self.members), len(self.groups)) 
	def __repr__(self):
		return '%s leaflet of a Membrane System with %s members and %s lipid groups' % (self.name, len(self.members), len(self.groups)) 

	def AddMember(self, index, resname):
		
		if	len(self.members) == 0:
			self.members.append([index, resname])
			self.groups.append(LipidGroup(resname))
			self.groups[0].AddMember(index)
			self.group_dict.update({resname:0})
		else:
			self.members.append([index, resname])
			addgroup = True
			group_ind = 0
			for rn in self.groups:
				if	resname == rn.lg_name:
					addgroup = False
					break
				group_ind+=1
			if	addgroup:
				self.groups.append(LipidGroup(resname))
				ng = len(self.groups)
				self.groups[ng-1].AddMember(index)
				self.group_dict.update({resname:ng-1})
			else:
				self.groups[group_ind].AddMember(index)
			
			#self.members=sorted(self.members,key=lambda self.members:self.members[1])
			
		return

	def GetGroupIndices(self, group_name):
		indices = []
		if group_name == "all":
			for ele in self.group_dict:
				gindex = self.group_dict[ele]
				indices += self.groups[gindex].lg_members
		elif group_name in self.group_dict:
			gindex = self.group_dict[group_name]
			indices = self.groups[gindex].lg_members
		else:
			#unkwown group name- print warning and use the default "all"
			print "!! Warning - request for unknown Lipid Group \'",group_name,"\' from the ",self.name," leaflet"
			print "!! using the default \"all\""
			for ele in self.group_dict:
				gindex = self.group_dict[ele]
				indices += self.groups[gindex].lg_members

		return list(indices)

	def GetMemberIndices(self):
		indices = []
		for ele in self.members:
			indices.append(ele[0])
		return list(indices)

	def HasGroup(self, group_name):
		return [group_name in self.group_dict]
		
class LipidGroup:
	def __init__(self, name):
		self.lg_members = []
		self.lg_name = name
		return

	def AddMember(self, new_mem):
		self.lg_members.append(new_mem)
		return

## this is the main class - the Membrane System (MemSys) object
class MemSys:
	# pass the mda anaylis trajectory object and a selection with the membrane (i.e. w/o water and ions)
	# optional - specify the plane that the membrane is in - default is xy with normal in z
	def __init__(self, mda_traj, mem_sel, plane="xy",fskip=1):
		#defaults - xy plane with z normal
		ii=0
		jj=1
		kk=2	
		if	plane=="yz" or plane=="zy":
			ii=1
			jj=2
			kk=0
		if	plane=="xz" or plane=="zx":
			ii=0
			jj=2
			kk=1
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
		self.frame = frames
		#loop over the frames
		f=0
		for frame in mda_traj[::fskip]:
			#add the frame object for this frame
			self.frame.append(Frame(self.nlipids))
			# set the box dimensions and the time for this frame
			self.frame[f].SetBox(frame.dimensions[0:3])
			self.frame[f].SetTime(frame.time)
			# loop over the residues (lipids) and get the centers of mass
			r=0			
			for res in mem_sel.residues:
				self.frame[f].lipidcom[r].extract(res)
				r+=1
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
			#print "unwrapping frame ",frame.frame
			currcoord = frame._pos[index]
			if firstframe:
				oldcoord = np.copy(currcoord)
				firstframe = False
			else:
				abc = frame.dimensions[0:3]
				wrapcoord = mda_wrap_coordinates(abc, currcoord, oldcoord)
				frame._pos[index] = wrapcoord[:]
				oldcoord = np.copy(wrapcoord)
			r=0			
			for res in mem_sel.residues:
				self.frame[f].lipidcom[r].extract(res, unwrap=True)
				r+=1
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
		f=0
		for fr in self.frame:
			count=0
			for i in indices:
				com_curr = fr.lipidcom[i].com_unwrap[self.plane]
				selcoords[f,count]=com_curr[:]
				count+=1
			f+=1
		
		#initialize a numpy array to hold the msd for the selection		
		msd = np.zeros((self.nframes, 6))
		#initialize a running stats object to do the averaging
		drs_stat = RunningStats()
		#initialize a running stats object for the diffusion constant (frame/time average)
		diff_stat = RunningStats()
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
			for	val in drs:
				drs_curr = val[:]	
				drs_mag = drs_curr.sum()
				drs_stat.Push(drs_mag)
			#get the msd for the current selection
			msdcurr = drs_stat.Mean()
			devcurr = drs_stat.Deviation()
			#dt = times[i]-times[0]
			DiffCon = msdcurr/(2.0*dim*dt)
			diff_stat.Push(DiffCon)
			#print "msdcurr ",msdcurr
			#push to the msd array
			
			msd[i,0]=dt
			msd[i,1]=msdcurr
			msd[i,2]=devcurr
			msd[i,3]=DiffCon
			msd[i,4]=diff_stat.Mean()
			msd[i,5]=diff_stat.Deviation()
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
		
		for fr in self.frame:
			n=0
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
					if	rxy<distxy:
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
					if	rxy<distxy:
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
			f+=1
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
		xyzout=False
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
			ctime = self.frame[f].time
			clusters = []
#			masterlistf = []
#			masterlistf += masterlist
			#rebuild the master list each frame
			masterlistf = list()
			for i in indices:
				masterlistf.append([i, False])
#			print "master ",masterlistf
			boxc=self.frame[f].box
			boxc_xh = boxc[xi]/2.0
			boxc_yh = boxc[yi]/2.0
			#print boxc
			clustind = 0
			neighborlist = []
			while len(masterlistf)>0:
				#print "master ",masterlistf
				start = masterlistf[0][0]
				masterlistf[0][1]=True
			#	print 
				# reset the neighborlist
				neighborlist = []
				#seed the neighborlist with the start
				neighborlist.append(start)
				#now loop over the neighborlist and build neighbors and neighbors of neigbors for this cluster
				i=0
				while i < len(neighborlist):
					ele = neighborlist[i]
					startn = ele
					coms = self.frame[f].lipidcom[startn].com		
					#get neighbors of the start
					#mindex=0
					for j in xrange(len(masterlistf)):
					#for elem in masterlistf:
						elem = masterlistf[j]
						incluster = elem[1]
					#	print "second incluster ",incluster
						if not incluster:
							ci = elem[0]
							comc = self.frame[f].lipidcom[ci].com
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
							if	rxy <= dist:
								#print "passed! adding ",masterlistf[mindex][0]," to the neighborlist"
								neighborlist.append(masterlistf[j][0])
								masterlistf[j][1]=True
						#mindex+=1
					i+=1
				#filter the masterlistf
			#	print "neighlist", neighborlist
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
				if	size>maxi:
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
#			if	f == 0:
#				break
			if xyzout:
				# Open up the file to write to
				xyz_name = "clusters_frame"+str(f)+".xyz"
				xyz_out = open(xyz_name, "w")
				nats=0
				for cluster in clusters:
					nats+=len(cluster)
				xyz_out.write(str(nats))
				xyz_out.write("\n")
				xyz_out.write("clusters")
				xyz_out.write("\n")
				c=0
				for cluster in clusters:
					for index in cluster:
							line = str(c)+" "+str(self.frame[f].lipidcom[index].com[0])+" "+str(self.frame[f].lipidcom[index].com[1])+" "+str(self.frame[f].lipidcom[index].com[2])
							xyz_out.write(line)
							xyz_out.write("\n")	
					c+=1
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

	# function to compute the mean squared displace (msd) along with the diffusion constant of a group 
	def CalcAreaPerLipid(self, leaflet="both",group="all"):
		
		#diffusion dimension - assume lateral so, dim=2
		dim=2
		do_leaflet = []
		if leaflet == "both":
			do_leaflet.append('upper')
			do_leaflet.append('lower')
		
		elif leaflet == "upper" or leaflet == "lower":
			do_leaflet.append(leaflet)
		else:
			#unknown option--use default "both"
			print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the ",self.name," leaflet"
			print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
			
		xi = self.plane[0]
		yi = self.plane[1]
		zi = self.norm
		sub_fact = (2.0*np.pi/3.0 - np.sqrt(3.0)/2.0)
		#initialize a numpy array to hold the msd for the selection		
		areas = np.zeros((self.nframes, 4))
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
			
		f=0
		#loop over the frames
		for fr in self.frame:
			dt = fr.time
			boxc=self.frame[f].box
			boxc_xh = boxc[xi]/2.0
			boxc_yh = boxc[yi]/2.0
			area_stat_config = RunningStats()
			#loop over the leaflets
			for leaflets in do_leaflet:
				indices = indices_leaflet[leaflets]
				all_mem = all_mem_leaflet[leaflets]
				#loop over the group indices in this leaflet
				for	index in indices:
					comc = self.frame[f].lipidcom[index].com[:]
					rdist_min = 10000.0
					#loop over the COMs of non group 
					#get all the leaflet members
					
					for a in all_mem:
						#print "a ",a
						if a != index:
							comn = self.frame[f].lipidcom[a].com[:]
							dx = comc[xi]-comn[xi]
							dy = comc[yi]-comn[yi]
						
							#Minimum image -- coordinates must be pre-wrapped 
							if np.absolute(dx) > boxc_xh:
								dx = boxc[xi] - np.absolute(comc[xi]-boxc_xh) - np.absolute(comn[xi]-boxc_xh)
							if np.absolute(dy) > boxc_yh:
								dy = boxc[yi] - np.absolute(comc[yi]-boxc_yh) - np.absolute(comn[yi]-boxc_yh)
							rxy = np.sqrt(dx*dx+dy*dy)
							#print "rxy ",rxy," dx ",dx," dy ",dy
							if	rxy < rdist_min:
								rdist_min = rxy
					#got the min dist, now compute area
					#print "rdist_min ",rdist_min
					area = np.pi*rdist_min*rdist_min - (rdist_min*rdist_min)*sub_fact
					area_stat_config.Push(area)
			area_conf_avg = area_stat_config.Mean()
			area_stat.Push(area_conf_avg)
			area_time_run = area_stat.Mean()
			area_time_run_dev = area_stat.Deviation()
			areas[f][0]=dt
			areas[f][1]=area_conf_avg
			areas[f][2]=area_time_run
			areas[f][3]=area_time_run_dev
			f+=1
		return areas
		
		

	#Outputs the lipid com coordinates as a xyz trajectory file
	def WriteXYZ(self, xyz_name, unwrap=False):
		# Open up the file to write to
		xyz_out = open(xyz_name, "w")
		
		for f in self.frame:
			xyz_out.write(str(self.nlipids))
			xyz_out.write("\n")
			xyz_out.write("Membrane System")
			xyz_out.write("\n")
			l = 0
			for com in f.lipidcom:
				#write to file
				
				if	unwrap:
					lc = "L"
					if	bool(int(self.leaflet[l])):
						lc = "U"
					line = lc+str(com.type[0])+" "+str(com.com_unwrap[0])+" "+str(com.com_unwrap[1])+" "+str(com.com_unwrap[2])
				else:
					line = str(com.type[0])+" "+str(com.com[0])+" "+str(com.com[1])+" "+str(com.com[2])
				xyz_out.write(line)
				xyz_out.write("\n")
				l+=1
		xyz_out.close()


#	def FilterList(list_in):
#		l = len(list_in)-1
#		for i, v in enumerate(reversed(list_in)):
#			if v[1] == True:
#				list_in.pop(l-i)
#		return list_in

#	def FilterList(list_in):
#		
#		return [v for v in list_in if !v[1]]
