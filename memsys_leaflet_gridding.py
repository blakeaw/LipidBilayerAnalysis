'''
	Classes and functions to implement lipid COM gridding and analysis for lipid bilayers. Acts on MemSys objects.
	The gridding and anlaysis procedures are based on 
	the decriptions given in Gapsys et al. J Comput Aided Mol Des (2013) 27:845-858,
	which is itself a modified version of the GridMAT-MD method by Allen et al. Vol. 30, No. 12 Journal of Computational Chemistry.
	However, I have currently left out several bits of the extra functionality, e.g. the handling of an embedded protein. 
'''

import numpy as np
#import my running stats class
from RunningStats import *

class LipidGrid_2d:
	def __init__(self, ms_frame, ms_com_indices,ms_plane,nxbins=50,nybins=50):
		#store the frame and leaflet
		self.frame = ms_frame
		#self.leaflet = ms_leaflet
		#get the x and y indices
		ix = ms_plane[0]
		iy = ms_plane[1]
		iz = [i for i in [0,1,2] if i not in ms_plane][0]
		#get the box dimemsions
		box = ms_frame.box[ms_plane]
		boxx = box[ix]
		boxy = box[iy]
		#save the numbers of bins
		self.x_nbins = nxbins
		self.y_nbins = nybins
		#initialize the edges of the and centers of the gridpoints
		# x
		self.x_min = 0.0
		self.x_max = boxx
 		self.x_edges = np.linspace(self.x_min,self.x_max,(nxbins+1),endpoint=True)
		self.x_incr = self.x_edges[1]-self.x_edges[0]
		x_incr_h = self.x_incr/2.0
		self.x_centers = np.zeros(nxbins)
		self.x_nedges = len(self.x_edges)
		for i in xrange(1,self.x_nedges):
			j=i-1
			self.x_centers[j]=self.x_edges[j]+x_incr_h

		# y
		self.y_min = 0.0
		self.y_max = boxy
 		self.y_edges = np.linspace(self.y_min,self.y_max,(nybins+1),endpoint=True)
		self.y_incr = self.y_edges[1]-self.y_edges[0]
		y_incr_h = self.y_incr/2.0
		self.y_centers = np.zeros(nybins)
		self.y_nedges = len(self.y_edges)
		for i in xrange(1,self.x_nedges):
			j=i-1
			self.y_centers[j]=self.y_edges[j]+y_incr_h
		self.x_length = self.x_max-self.x_min
		self.y_length = self.y_max-self.y_min
		# get the lipid indices for this leaflet
		indices = ms_com_indices
		#now assign lipids to the gridpoints
		self.lipid_grid = np.zeros((nxbins,nybins),dtype=np.int)
		self.lipid_grid_z = np.zeros((nxbins,nybins))
		bxh = boxx/2.0
		byh = boxy/2.0
		cx = 0
		for x in self.x_centers:
			cy = 0
			for y in self.y_centers:
				r_min = 1.0e10
				i_min = 0
				z_min = 0.0
				for i in indices:
					xi = ms_frame.lipidcom[i].com[ix]
					yi = ms_frame.lipidcom[i].com[iy]
					zi = ms_frame.lipidcom[i].com[iz]
					#print "iz ",iz," zi ",zi
					dx = x-xi
					dy = y-yi
					#Minimum image -- coordinates must be pre-wrapped 
					if np.absolute(dx) > bxh:
						dx = boxx - np.absolute(x-bxh) - np.absolute(xi-bxh)
					if np.absolute(dy) > bxh:
						dy = boxy - np.absolute(y-byh) - np.absolute(yi-byh)
					rxy = np.sqrt(dx**2+dy**2)
					if rxy < r_min:
						r_min=rxy
						i_min = i
						z_min = zi
				#print "i_min ",i_min," z_min ",z_min
				self.lipid_grid[cx,cy]=i_min
				self.lipid_grid_z[cx,cy]=z_min
				cy+=1
			cx+=1
	
	def GetIndexAt(self,ix,iy):
		return self.lipid_grid[ix,iy]

	def GetZAt(self,ix,iy):
		return self.lipid_grid_z[ix,iy]
			
	#Outputs the grid as an xyz coordinate file
	def WriteXYZ(self, xyz_name):
		# Open up the file to write to
		xyz_out = open(xyz_name, "w")
		npoints = self.x_nbins*self.y_nbins
		comment = "Leaflet Grid "+self.leaflet.name
		xyz_out.write(str(npoints))
		xyz_out.write("\n")
		xyz_out.write(comment)
		xyz_out.write("\n")
		
		cx=0
		for x in self.x_centers:
			cy=0
			for y in self.y_centers:
				#get the z coordinate
				z = self.lipid_grid_z[cx,cy]
				#get the lipid resname
				ic = self.lipid_grid[cx,cy]
				oname = self.frame.lipidcom[ic].type
				#write to file
				
			
				line = str(oname)+" "+str(x)+" "+str(y)+" "+str(z)
			
				xyz_out.write(line)
				xyz_out.write("\n")
				cy+=1
			cx+=1			
		xyz_out.close()
		return

class LeafletGrids:
	def __init__(self, ms_frame, ms_leaflets,ms_plane,nxbins=50,nybins=50):
		#store the frame and leaflet
		self.frame = ms_frame
		self.com_leaflets = ms_leaflets
		self.plane = ms_plane
		self.norm = [i for i in [0,1,2] if i not in ms_plane][0]
		self.nbins_x = nxbins
		self.nbins_y = nybins
		self.leaflets = {}
		self.myframe = ms_frame.number
		#initialize the grids
		#upper
		upper_indices = ms_leaflets['upper'].GetMemberIndices()
		self.leaflets['upper'] = LipidGrid_2d(ms_frame,upper_indices,ms_plane,nxbins=nxbins,nybins=nybins)
		#lower
		lower_indices = ms_leaflets['lower'].GetMemberIndices()
		self.leaflets['lower'] = LipidGrid_2d(ms_frame,lower_indices,ms_plane,nxbins=nxbins,nybins=nybins)
		return

	def ThicknessGrid(self):
		tgrid = np.zeros((self.nbins_x,self.nbins_y))
		for ix in xrange(self.nbins_x):
			for iy in xrange(self.nbins_y):
				zu = self.leaflets['upper'].GetZAt(ix,iy) 
				zl = self.leaflets['lower'].GetZAt(ix,iy)
				dz = zu - zl
				tgrid[ix,iy]=dz 
		return tgrid
	
	def AverageThickness(self,return_grid=False):
		trun = RunningStats()
		tgrid = self.ThicknessGrid()
		for ix in xrange(self.nbins_x):
			for iy in xrange(self.nbins_y):
				tc = tgrid[ix,iy]
				trun.Push(tc) 
		avg_out = (trun.Mean(),trun.Deviation())
		if return_grid:
			return avg_out,tgrid
		else:
			return avg_out

	def MapToGrid(self,com_values_dict,leaflet='both'):
		do_leaflet = []
		if leaflet == "both":
			do_leaflet.append('upper')
			do_leaflet.append('lower')
		elif leaflet == "upper" or leaflet == "lower":
			do_leaflet.append(leaflet)
		else:
			#unknown option--use default "both"
			print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the LeafletGrids of frame ",self.myframe
			print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
			do_leaflet.append('upper')
			do_leaflet.append('lower')
		
		out_dict = {}		
		for leaf in do_leaflet:
			out_dict[leaf] = np.zeros((self.nbins_x,self.nbins_y))
			for ix in xrange(self.nbins_x):
				for iy in xrange(self.nbins_y):
					com_ind=self.leaflets[leaf].GetIndexAt(ix,iy)
					value = com_values_dict[com_ind]
					out_dict[leaf][ix,iy]=value
					
		return out_dict
	
	def AreaPerLipid(self):
		do_leaflet = []
		do_leaflet.append('upper')
		do_leaflet.append('lower')
		#get the unique type/resnames in the system
		resnames = []
		for leaf in do_leaflet:
			for group in self.com_leaflets[leaf].groups:
				gname = group.name()
				if gname not in resnames:
					resnames.append(gname)
		#initialize counters for each residue/type
		area_run_per_res_type = {}
		
		for name in resnames:
			area_run_per_res_type[name]=RunningStats()
			
		area_per_lipid = {}
		
		area_run = RunningStats()
		for leaf in do_leaflet:
			area_per_bin = self.leaflets[leaf].x_incr*self.leaflets[leaf].y_incr
			lip_ind = self.com_leaflets[leaf].GetMemberIndices()
			for i in lip_ind:
				rname = self.frame.lipidcom[i].type
				locations = np.where(self.leaflets[leaf].lipid_grid == i)
				nlocs = len(locations[0])
				#print locations
				#print 'nlocs ',nlocs
				area = area_per_bin*nlocs
				area_per_lipid[i]=area
				area_run_per_res_type[rname].Push(area)
				
				area_run.Push(area)
		
		average_per_res = {}
		for name in resnames:
			average = area_run_per_res_type[name].Mean()	
			std = area_run_per_res_type[name].Deviation()
			average_per_res[name] = (average,std)		
		system_average = area_run.Mean()
		system_dev = area_run.Deviation()			
		
		output = ((system_average,system_dev),average_per_res,area_per_lipid)		
		return output
#	def Curvature(self):
		

	def GridToDict(self,in_grid,leaflet='upper'):
		out_dict = {}
		for ix in xrange(self.nbins_x):
			for iy in xrange(self.nbins_y):
				l_i = self.leaflets[leaflet].GetIndexAt(ix,iy) 
				grid_val = in_grid[ix,iy]
				out_dict[l_i]=grid_val		 
		return out_dict

	def GetXYZ(self,leaflet='both',zvalue_dict='Default',color_dict=None):
		do_leaflet = []
		if leaflet == "both":
			do_leaflet.append('upper')
			do_leaflet.append('lower')
		elif leaflet == "upper" or leaflet == "lower":
			do_leaflet.append(leaflet)
		else:
			#unknown option--use default "both"
			print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the LeafletGrids of frame ",self.myframe
			print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
			do_leaflet.append('upper')
			do_leaflet.append('lower')
		out_dict = {}
		npoints = (self.nbins_x*self.nbins_y)
		X = np.zeros(npoints)
		Y = np.zeros(npoints)
		Z = np.zeros(npoints)
		C = np.zeros(npoints)
		for leaf in do_leaflet:
			npt = 0
			cx=0
			for x in self.leaflets[leaf].x_centers:
				cy=0
				for y in self.leaflets[leaf].y_centers:
					#get the z coordinate
					z = self.leaflets[leaf].GetZAt(cx,cy)
					ic = self.leaflets[leaf].GetIndexAt(cx,cy)
					#optionally pull z value from lipid index dictionary
					if zvalue_dict is not 'Default':
						z = zvalue_dict[ic] 
					X[npt]=x
					Y[npt]=y
					Z[npt]=z
					if color_dict is not None:
						C[npt]=color_dict[ic]
					npt+=1
					cy+=1
				cx+=1	
			if color_dict is not None:
				col_min = min(C)
				C-=col_min
				col_max = max(C)
				C/=col_max
				out_dict[leaf]=(X,Y,Z,C)
			else:	
				out_dict[leaf]=(X,Y,Z)
		return	out_dict
	
	def WriteXYZ(self,leaflet='both',zvalue_dict='Default',out_path="./"):
		do_leaflet = []
		if leaflet == "both":
			do_leaflet.append('upper')
			do_leaflet.append('lower')
		elif leaflet == "upper" or leaflet == "lower":
			do_leaflet.append(leaflet)
		else:
			#unknown option--use default "both"
			print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the LeafletGrids of frame ",self.myframe
			print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
			do_leaflet.append('upper')
			do_leaflet.append('lower')
		out_name = out_path+"leaflet_grid_f"+str(self.myframe)+"_"
		for leaf in do_leaflet:
			out_name+=leaf[0]
		out_name+=".xyz"	
		# Open up the file to write to
		xyz_out = open(out_name, "w")
		npoints = (self.nbins_x*self.nbins_y)*len(do_leaflet)
		comment = "Leaflet Grid in xyz coordinate format"
		xyz_out.write(str(npoints))
		xyz_out.write("\n")
		xyz_out.write(comment)
		xyz_out.write("\n")
		
		for leaf in do_leaflet:
			cx=0
			for x in self.leaflets[leaf].x_centers:
				cy=0
				for y in self.leaflets[leaf].y_centers:
					#get the z coordinate
					z = self.leaflets[leaf].GetZAt(cx,cy)
					ic = self.leaflets[leaf].GetIndexAt(cx,cy)
					#optionally pull z value from lipid index dictionary
					if zvalue_dict is not 'Default':
						z = zvalue_dict[ic] 
					#get the lipid resname
					
					oname = self.frame.lipidcom[ic].type
					#write to file
				
			
					line = str(oname)+" "+str(x)+" "+str(y)+" "+str(z)
			
					xyz_out.write(line)
					xyz_out.write("\n")
					cy+=1
				cx+=1			
		xyz_out.close()
		return




