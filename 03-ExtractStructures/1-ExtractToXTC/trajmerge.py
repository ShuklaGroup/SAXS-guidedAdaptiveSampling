#!/home/czhao37/anaconda3/bin/python

import mdtraj as md
import argparse

def parse_cmdln():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input', dest='trajList',
                        help='File containing trajectory names.', required=True)
    parser.add_argument('-t', '--topology', dest='top',
                        help='File containing topology.', default=None)
    parser.add_argument('-o', '--output', dest='out',
                        help='Name of output trajectory.', default='trajout.xtc')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    options = parse_cmdln()
    file = open(options.trajList,"r")
    trajnames = [ line.rstrip() for line in file ]
    trajnum   = len(trajnames)
    print("Starting merging the trajectories: " + str(trajnum) + " in total ..." )
    trajfinal = md.load(trajnames[0], top=options.top)
    if trajnum > 0:
        trajnames = trajnames[1:]
        for trajname in trajnames:
            try:
                traj = md.load(trajname, top=options.top)
                trajfinal = trajfinal.join(traj)
            except:
                continue
        print("Starting saving the merged trajectories ...")
        trajfinal.save(options.out)

