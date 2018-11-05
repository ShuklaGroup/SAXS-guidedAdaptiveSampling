### For each cluster, extracting 100 structures randomly from the raw MD trajectories 
### to calculate the average SAXS profile of each state and the SAXS discrepancy scores
### between the target and each state.  
### Required packages: pickle, glob, msmbuilder, mdtraj
### Usage: python extractClusterToPDB.py (-h for options)
### @Chuankai Zhao, czhao37@illinois.edu

from msmbuilder.utils import io
import pickle
import glob
import mdtraj as md
import argparse

### Given the frame indexes, load and save as pdb for the SAXS calculation. 
def pickstates(selects, trajnames, top, stride, NP):
    select   = selects[options.NP]
    cl_state = NP
    count = 0
    for structure in select:
        filename = "state" + str(cl_state) + '_' + str(count) + '.pdb'
        traj     = md.load(trajnames[structure[0]], top=top, frame=structure[1]*stride)
        traj.save_pdb(filename)
        count    = count + 1

### Define command line options
def parse_cmdln():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-m', '--msm-file', dest='msm',
                        help='MSM File.', default='MSM.pkl', required=True)
    parser.add_argument('-c', '--cluster-file', dest='cl',
                        help='Cluster File.',
                        default='clustering.pkl', required=True)
    parser.add_argument('-ns', '--n-samples', dest='NS',
                        help='Number of frames to extract.',
                        default=100, type=int)
    parser.add_argument('-np', '--n-processors', dest='NP',
                        help='State ID.',
                        default=10, type=int)
    parser.add_argument('-s', '--stride', dest='stride',
                        help='Stride of the raw trajectories being subsampled.', default=None, type=int, required=True)
    args = parser.parse_args()
    return args

### Main program
if __name__ == '__main__':
    options   = parse_cmdln()
    msm_name  = options.msm
    cl_name   = options.cl
    cluster   = pickle.load(open(cl_name,'rb'), encoding='latin1')
    msm       = io.load(msm_name)
    cl_label  = cluster.labels_
    n_samples = options.NS
    stride    = options.stride

    # Read the topology file and the list of MD raw trajectories
    top       = "/home/amoffet2/msm_network_project/folding/lindorff-larsen_2011_trajs/protein_g-350K/protein_g.pdb"
    trajnames = glob.glob("/home/amoffet2/msm_network_project/folding/lindorff-larsen_2011_trajs/protein_g-350K/DESRES-Trajectory_NuG2-*-protein/NuG2-*-protein/*.dcd")

    # Using msmbuilder.draw_samples function to extract 100 structures from each state
    selects   = msm.draw_samples(cl_label,n_samples)
    pickstates(selects, trajnames, top, stride, options.NP)
