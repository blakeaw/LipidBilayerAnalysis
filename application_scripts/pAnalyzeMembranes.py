import MDAnalysis as mda
import numpy as np
import MemSys_v2 as ms
#import matplotlib.pyplot as plt
import pickle

psf='step5_charmm2amber.psf'
traj_in_pre="step7_"
traj_in_post=".nc"
traj_start = 1
traj_end = 5
frame_skip = 5
group_list = ['DOPE','POPC']
#group_list = ['DOPE','POPC','TLCL2']
outputpath='/home/wi8/LipidOut/Analysis_Outputs/mem0/'


print "Creating mda universe"
u = mda.Universe(psf)

membrane = u.select_atoms("not resname CLA and not resname TIP3 and not resname POT")

traj_list = []
print "looping over trajectories"
for	i in xrange(traj_start,(traj_end+1)):
	traj_name=traj_in_pre+str(i)+traj_in_post
	print "Doing trajectory ",traj_name
	traj_list.append(traj_name)

#load the trajectories
u.load_new(traj_list)

print "building membrane system object"
mem_sys = ms.MemSys(u.trajectory,membrane,fskip=frame_skip,frame_path=outputpath)
#name = outputpath+"mem_sys_bak.pickle"
#outfile = open(name,'wb')
#pickle.dump(mem_sys,outfile)
#outfile.close()
print "measuring MSD..."

msd_all = mem_sys.CalcMSD()
name = outputpath+"msd_all.pickle"
outfile = open(name,'wb')
pickle.dump(msd_all,outfile)
outfile.close()
#np.savetxt('msd_all.dat',msd_all)
for item in group_list:
	msd = mem_sys.CalcMSD(group=item)
	msd_upper = mem_sys.CalcMSD(leaflet='upper',group=item)
	msd_lower = mem_sys.CalcMSD(leaflet='lower',group=item)
	name1 = outputpath+"msd_"+item+".pickle"
	name2 = outputpath+"msd_"+item+"_upper.pickle"
	name3 = outputpath+"msd_"+item+"_lower.pickle"
	outfile = open(name1,'wb')
	pickle.dump(msd,outfile)
	outfile.close()
	outfile = open(name2,'wb')
	pickle.dump(msd_upper,outfile)
	outfile.close()
	outfile = open(name3,'wb')
	pickle.dump(msd_lower,outfile)
	outfile.close()
#	np.savetxt(name1,msd)
#	np.savetxt(name2,msd_upper)
#	np.savetxt(name3,msd_lower)

print "measuring membrane thickness..."
thickness,membrane_maps = mem_sys.CalcMembraneThickness()
#np.savetxt('mem_thickness.dat',thickness)
#np.savetxt('mem_thickness_map.dat',membrane_maps)
name = outputpath+"mem_thickness.pickle"
outfile = open(name,'wb')
pickle.dump(thickness,outfile)
outfile.close()
name = outputpath+"mem_thickness_map.pickle"
outfile = open(name,'wb')
pickle.dump(membrane_maps,outfile)
outfile.close()
print "measuring areas per lipid..."
area_all = mem_sys.CalcAreaPerLipid()
name = outputpath+"apl_all.pickle"
outfile = open(name,'wb')
pickle.dump(area_all,outfile)
outfile.close()
#np.savetxt('apl_all.dat',area_all)
for item in group_list:
	apl = mem_sys.CalcAreaPerLipid(group=item)
	apl_upper = mem_sys.CalcAreaPerLipid(leaflet='upper',group=item)
	apl_lower = mem_sys.CalcAreaPerLipid(leaflet='lower',group=item)
	name1 = outputpath+"apl_"+item+".pickle"
	name2 = outputpath+"apl_"+item+"_upper.pickle"
	name3 = outputpath+"apl_"+item+"_lower.pickle"
	outfile = open(name1,'wb')
	pickle.dump(apl,outfile)
	outfile.close()
	outfile = open(name2,'wb')
	pickle.dump(apl_upper,outfile)
	outfile.close()
	outfile = open(name3,'wb')
	pickle.dump(apl_lower,outfile)
	outfile.close()
	#np.savetxt(name1,apl)
	#np.savetxt(name2,apl_upper)
	#np.savetxt(name3,apl_lower)
#check clustering for cardiolipin - if it is in the system
if 'TLCL2' in group_list:
	print "Checking cardiolipin clusters..."
	cluster_avgs_upper = mem_sys.CheckClustering(leaflet='upper',group='TLCL2')
	cluster_map_upper = mem_sys.ExportClustersForPlotting()
	cluster_avgs_lower = mem_sys.CheckClustering(leaflet='lower',group='TLCL2')
	cluster_map_lower = mem_sys.ExportClustersForPlotting()
	cluster_avgs_both = mem_sys.CheckClustering(group='TLCL2')
	cluster_map_both = mem_sys.ExportClustersForPlotting()
	name1 = outputpath+"cardio-clust_upper.pickle"
	name2 = outputpath+"cardio-clust_lower.pickle"
	name3 = outputpath+"cardio-clust_both.pickle"
	outfile = open(name1,'wb')
	pickle.dump(cluster_avgs_upper,outfile)
	outfile.close()
	outfile = open(name2,'wb')
	pickle.dump(cluster_avgs_lower,outfile)
	outfile.close()
	outfile = open(name3,'wb')
	pickle.dump(cluster_avgs_both,outfile)
	outfile.close()
	name1 = outputpath+"cardio-clust_map-upper.pickle"
	name2 = outputpath+"cardio-clust_map-lower.pickle"
	name3 = outputpath+"cardio-clust_map-both.pickle"
	outfile = open(name1,'wb')
	pickle.dump(cluster_map_upper,outfile)
	outfile.close()
	outfile = open(name2,'wb')
	pickle.dump(cluster_map_lower,outfile)
	outfile.close()
	outfile = open(name3,'wb')
	pickle.dump(cluster_map_both,outfile)
	outfile.close()
#plt.plot(msd_dope[:,0],msd_dope[:,4])
#plt.show()
