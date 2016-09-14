import MDAnalysis as mda
import numpy as np
import LipidBilayerAnalysis.MemSys as ms
import LipidBilayerAnalysis.memsys_leaflet_gridding as lg
import LipidBilayerAnalysis.plot_generation_functions as pgf
import LipidBilayerAnalysis.mda_density_profile as dp
import LipidBilayerAnalysis.mda_deuterium_order_parameter as dop
#import matplotlib.pyplot as plt
import pickle
from matplotlib import pyplot as plt

def dot_to_d(string):
    new_string = ''
    for c in string:
        if c == '.':
            new_string+='d'
        else:
            new_string+=c
    return new_string


psf = 'step5_charmm2amber.psf'
traj_in_pre = "step7_"
traj_in_post = ".nc"
traj_start = 1
traj_end = 3
ns_per_traj = 10.0
frame_per_traj = 500
frame_skip = 10
group_list = ['DOPE','POPC']
#group_list = ['DOPE','POPC','TLCL2']
outputpath = '/home/wi8/LipidOut/Analysis_Outputs/mem0/'
nprocs = 10

print "Creating mda universe"
u = mda.Universe(psf)

bilayer = u.select_atoms("not resname CLA and not resname TIP3 and not resname POT")

traj_list = []
print "looping over trajectories"
for    i in xrange(traj_start,(traj_end+1)):
    traj_name=traj_in_pre+str(i)+traj_in_post
    print "Doing trajectory ",traj_name
    traj_list.append(traj_name)

#load the trajectories
u.load_new(traj_list)

print "building bilayer system object"
mem_sys = ms.MemSys(u.trajectory,bilayer,fskip=frame_skip,frame_path=outputpath,nprocs=4)
#remove leaflet COM motion
print "removing the leaflet COM motion"
mem_sys.RemoveLeafletCOMmotion()


print "measuring MSD..."

msd_dat = []
msd_names = []
msd_all = mem_sys.CalcMSD_parallel(nprocs=nprocs)
name = outputpath+"msd_all.pickle"
outfile = open(name,'wb')
pickle.dump(msd_all,outfile)
outfile.close()
msd_dat.append(msd_all.copy())
msd_names.append('All')
del msd_all

for item in group_list:
    msd = mem_sys.CalcMSD_parallel(group=item,nprocs=nprocs)
    name1 = outputpath+"msd_"+item+".pickle"
    outfile = open(name1,'wb')
    pickle.dump(msd,outfile)
    outfile.close()
    msd_dat.append(msd.copy()) 
    if item == 'TLCL2':
        msd_names.append('CL')
    else:
        msd_names.append(item)

pgf.plot_msd(msd_dat,name_list=msd_names,filename=outputpath+'msd.eps')
pgf.plot_msd(msd_dat,name_list=msd_names,filename=outputpath+'msd.png')
del msd_dat
del msd_names

#print "measuring bilayer thickness..."
#thickness_out = mem_sys.CalcMembraneThickness_parallel(nprocs=nprocs)
#thickness = thickness_out[0]
#membrane_maps = thickness_out[1]
#name = outputpath+"mem_thickness.pickle"
#outfile = open(name,'wb')
#pickle.dump(thickness,outfile)
#outfile.close()
#name = outputpath+"mem_thickness_map.pickle"
#outfile = open(name,'wb')
#pickle.dump(membrane_maps,outfile)
#outfile.close()
#del thickness_out
#del thickness
#del membrane_maps

print "measuring areas per lipid..."
area_all = mem_sys.CalcAreaPerLipid_Box()
name = outputpath+"apl_all.pickle"
outfile = open(name,'wb')
pickle.dump(area_all,outfile)
outfile.close()

apl_list = [ area_all.copy() ]
apl_names = [ 'All' ]
pgf.plot_area_per_lipid(apl_list,filename=outputpath+'apl_box.eps',interval=10)
pgf.plot_area_per_lipid(apl_list,filename=outputpath+'apl_box.png', interval=10)
for item in group_list:
    apl = mem_sys.CalcAreaPerLipid_ClosestNeighborCircle(group=item)
    name1 = outputpath+"apl_"+item+".pickle"
    outfile = open(name1,'wb')
    pickle.dump(apl,outfile)
    outfile.close()
    apl_list.append(apl.copy()) 
    if item == 'TLCL2':
        apl_names.append('CL')
    else:
        apl_names.append(item)

pgf.plot_area_per_lipid(apl_list,name_list=apl_names, interval=10)
pgf.plot_area_per_lipid(apl_list,name_list=apl_names,filename=outputpath+'apl.png', interval=10)
del apl_list
del apl_names
del area_all

#now do displacement vectors
print "computing displacement vectors..."

frstep = 10
disp_upper = mem_sys.StepVector(leaflet='upper', fstep=frstep, wrapped=True)
disp_u_colors = mem_sys.StepVectorColors(leaflet='upper')
outfile = open(outputpath+"dispvecs_upper.pickle",'wb')
pickle.dump(disp_upper, outfile)
outfile.close()
outfile = open(outputpath+"dispvecs_upper_colors.pickle",'wb')
pickle.dump(disp_u_colors, outfile)
outfile.close()

disp_lower = mem_sys.StepVector(leaflet='lower', fstep=frstep, wrapped=True)
disp_l_colors = mem_sys.StepVectorColors(leaflet='lower')
outfile = open(outputpath+"dispvecs_lower.pickle",'wb')
pickle.dump(disp_upper, outfile)
outfile.close()
outfile = open(outputpath+"dispvecs_lower_colors.pickle",'wb')
pickle.dump(disp_u_colors, outfile)
outfile.close()

ndu = len(disp_upper)
ndl = len(disp_lower)

ntraj = traj_end - traj_start +1
print "ntraj ",ntraj
total_ns = ntraj * ns_per_traj
print "total ns ",total_ns
total_frame = ntraj * frame_per_traj 
print "total_frame ",total_frame
ns_per_frame = total_ns / total_frame
print "ns per frame ",ns_per_frame
ns_per_ms_frame = ns_per_frame * frame_skip
print "ns per ms frame ", ns_per_ms_frame
ns_per_stepvec = ns_per_ms_frame*frstep
print "ns per stepvec ", ns_per_stepvec
# upper
for i in xrange( ndu ):
    j = i + 1
    t_one = ns_per_stepvec*i
    t_two = ns_per_stepvec*j
    t_one_s = '%.2f' % t_one
    t_two_s = '%.2f' % t_two
    t_one_s = dot_to_d(t_one_s)
    t_two_s = dot_to_d(t_two_s)
    name_eps = "dispvecs_upper_" + t_one_s + "t" + t_two_s + "ns.eps" 
    name_png = "dispvecs_upper_" + t_one_s + "t" + t_two_s + "ns.png" 
    vec_curr = disp_upper[i]
    pgf.plot_step_vectors(vec_curr, colors=disp_u_colors[0], filename=outputpath+name_eps)
    pgf.plot_step_vectors(vec_curr, colors=disp_u_colors[0], filename=outputpath+name_png)
# lower
for i in xrange( ndl ):
    j = i + 1
    t_one = ns_per_stepvec*i
    t_two = ns_per_stepvec*j
    t_one_s = '%.2f' % t_one
    t_two_s = '%.2f' % t_two
    t_one_s = dot_to_d(t_one_s)
    t_two_s = dot_to_d(t_two_s)
    name_eps = "dispvecs_lower_" + t_one_s + "t" + t_two_s + "ns.eps" 
    name_png = "dispvecs_lower_" + t_one_s + "t" + t_two_s + "ns.png" 
    vec_curr = disp_lower[i]
    pgf.plot_step_vectors(vec_curr, colors=disp_l_colors[0], filename=outputpath+name_eps)
    pgf.plot_step_vectors(vec_curr, colors=disp_l_colors[0], filename=outputpath+name_png)

del disp_lower
del disp_upper
del disp_l_colors
del disp_u_colors

#done with the COM representation calculations
del mem_sys

#electron density profiles
print "computing electron density profiles..."

sel_water = u.select_atoms('resname TIP3')
sel_watlip = u.select_atoms('not resname CLA and not resname POT')
sel_dope = u.select_atoms('resname DOPE')
sel_popc = u.select_atoms('resname POPC')
sel_cl = u.select_atoms('resname TLCL2')
fstep = 100
#dp_water = dp.ElectronDensityProfile_gaussians(u.trajectory, sel_water, fstep=fstep, refsel=bilayer, nbins=20)
dp_watlip = dp.ElectronDensityProfile_gaussians(u.trajectory, sel_watlip, fstep=fstep, refsel=bilayer, nbins=50)
#dp_bilayer = dp.ElectronDensityProfile_gaussians(u.trajectory, bilayer, fstep=fstep, refsel=bilayer, nbins=20)
#dp_dope = dp.ElectronDensityProfile_gaussians(u.trajectory, sel_dope, fstep=fstep, refsel=bilayer, nbins=20)
#dp_popc = dp.ElectronDensityProfile_gaussians(u.trajectory, sel_popc, fstep=fstep, refsel=bilayer, nbins=20)
plt.plot(dp_watlip[0], dp_watlip[1], linewidth=2.0, label='Water+Bilayer')
#plt.plot(dp_water[0], dp_water[1], linewidth=2.0, label='Water')
#plt.plot(dp_bilayer[0], dp_bilayer[1], linewidth=2.0, label='Bilayer')
#plt.plot(dp_dope[0], dp_dope[1], linewidth=2.0, label='DOPE')
#plt.plot(dp_popc[0], dp_popc[1], linewidth=2.0, label='POPC')
#if 'TLCL2' in group_list:
#    dp_cl = dp.ElectronDensityProfile_gaussians(u.trajectory, sel_cl, fstep=fstep, refsel=bilayer, nbins=20)
#    plt.plot(dp_cl[0], dp_cl[1], linewidth=2.0, label='CL')
plt.savefig(outputpath+'electron_density.eps')
plt.savefig(outputpath+'electron_density.png')

plt.show()

