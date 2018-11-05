### Featurization based on dihedral angles for the protein folding trajectories
### Required packages: mdtraj, msmbuilder, glob
### @Chuankai Zhao, czhao37@illinois.edu

import mdtraj as md
import glob
from msmbuilder.featurizer import DihedralFeaturizer
from msmbuilder.utils import verbosedump, verboseload

# Set the path of MD trajectories and the name of topology files.
trajaddress = "/home/amoffet2/msm_network_project/folding/lindorff-larsen_2011_trajs/protein_g-350K/DESRES-Trajectory_NuG2-*-protein/NuG2-*-protein/*.dcd"
top = "/home/amoffet2/msm_network_project/folding/lindorff-larsen_2011_trajs/protein_g-350K/protein_g.pdb"

# Load the trajectories using mdtraj
files = glob.glob(trajaddress)
traj_list=[]
for f in files:
   t = md.load(f, top=top, stride=10)
   traj_list.append(t)

# Featurize the trajectories based on phi, psi, chi1 dihedral angles
model = DihedralFeaturizer(types=['phi', 'psi','chi1'])
features = model.transform(traj_list)

# Set the path of output file and save the output. 
pkl = "/home/czhao37/2-SAXS-Adaptive_Samping/6-protein-G/features/featurized_1mio.pkl"
verbosedump(features, pkl)
