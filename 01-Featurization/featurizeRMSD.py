### Featurization based on RMSD from the native state for the protein folding trajectories
### Required packages: mdtraj, numpy
### @Chuankai Zhao, czhao37@illinois.edu

import mdtraj as md
import numpy as np

# Read the list of MD trajectories to featurize
trajnames = [ line.rstrip() for line in open("List","r") ]

# Read the native crystal structure
reference = md.load("native.pdb") 
top = reference.topology

# Calculate RMSDs of all frames in MD trajectories. 
for trajname in trajnames:
    traj = md.load( trajname + ".xtc", top=top, stride=8)
    rmsd = md.rmsd(traj, reference, atom_indices=top.select("name 'CA' or name 'N' or name 'O' or name 'CB' or name 'C'"))
    np.save(trajname + ".npy", rmsd)

