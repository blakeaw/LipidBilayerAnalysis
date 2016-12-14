# LipidBilayerAnalysis

Mean Squared displacement functions:
  pMSD.py - a set of functions to compute the MSD of MD Analysis selections and lists of selections
  pMSD_cyton.pyx - a Cython version of the MSD_list function from pMSD.py
  mda_msd.py - combines most of the features of the individual functions in pMSD.py into a single function

Center of Mass (COM) representation functions and classes:
  MemSys.py - processes the trajectory and selection of the membrane lipids and
	      rebuilds the trajectory with just the centers of mass of the lipids
	      Includes various anlyses based on the COM mass motion of the lipids.	
	Version 2 - The frames of COM representations are stored on disk using the shelve module.
	            This allows large trajectories to be processed without eating up all the RAM 
                    (although this requires up to a few GB of free disk space).
				
Dependencies:
    NumPy
    SciPy
    Matplotlib
    Seaborn
