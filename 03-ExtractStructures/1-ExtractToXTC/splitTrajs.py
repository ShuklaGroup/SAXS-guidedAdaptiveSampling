### Splitting individual trajectory into individual blocks according to the cluster index of each frame.
### Usage: python splitTrajs.py 
### Required packages: mdtraj, msmbuilder, numpy, pickle
### @Chuankai Zhao, czhao37@illinois.edu

import mdtraj as md
import pickle
import numpy as np

### load MD trajectory
def loadtraj(trajname):
    traj = md.load(trajname + ".mdcrd",top="ABA_Dimer.prmtop")
    return traj

### splitting individual trajectory into indidual blocks according to the cluster index of each frame.
def savetrajstate(trajname, traj, trajID):
    states  = cluster.labels_[trajID]
    frames  = np.shape(states)[0]
    stateID = states[0]
    frame = 0
    ID_serial = 0
    ID_Start = 0
    for state in states:
        if frame == frames - 1:
            if state == stateID:
                ID_End = frames
                ID_serial = ID_serial + 1
                tag = trajname + "_STATE_" + str(stateID) + "_" + str(ID_serial) + ".xtc"
                traj[ID_Start:ID_End].save_xtc(tag)
                ID_Start = frame

        if state != stateID:
            ID_End = frame
            ID_serial = ID_serial + 1
            tag = trajname + "_STATE_" + str(stateID) + "_" + str(ID_serial) + ".xtc"
            traj[ID_Start:ID_End].save_xtc(tag)
            stateID = state
            ID_Start = frame
            if frame == frames - 1:
                ID_End = frames
                ID_serial = ID_serial + 1
                tag = trajname + "_STATE_"  + str(stateID) + "_" + str(ID_serial) + ".xtc"
                traj[ID_Start:ID_End].save_xtc(tag)
                ID_Start = frame

        frame = frame + 1

### Main function
if __name__ == '__main__':

    ### Read the list of MD trajectories
    trajlist = "List"
    trajnames = [ line.rstrip() for line in open(trajlist, "r")]

    ### Read the cluster file
    cluster = pickle.load(open("clustering_tica.pkl","rb"))


    ### For each trajectory, splitting into individual blocks based on the cluster index
    trajID = 0
    for trajname in trajnames:
        traj = loadtraj(trajname)
        savetrajstate(trajname, traj, trajID)
        trajID = trajID + 1
