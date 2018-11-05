### Featurization based on the defined metrics (Zhao and Shukla, Sci. Rep., 2018) for the protein-protein association trajectories
### Features: dx, dy, dz, dAlpha, dBeta, dGamma, RMSD1, RMSD2
### Required packages: mdtraj, numpy
### @Chuankai Zhao, czhao37@illinois.edu

import mdtraj as md
import numpy as np

# Read the crystal structure for the reference structure in RMSD calculations
reference = md.load("ABA_Dimer.pdb")
top=reference.topology

# Read the list of MD trajectories
trajnames = [ line.rstrip() for line in open("List","r") ]

# Featurize the trajectories 
for trajname in trajnames:

        ### Load each trajectory
        traj = md.load(trajname + "_stripped.mdcrd", top="ABA_Dimer.prmtop")
        frames = traj.n_frames

        ### RMSDs of monomer 1 and monomer 2 (RMSD1, RMSD2) to consider the internal flexibility.
        reference_M1 = reference.atom_slice(top.select("resid 0 to 184"))
        reference_M2 = reference.atom_slice(top.select("resid 185 to 369"))
        rmsd_M1 = md.rmsd(traj.atom_slice(top.select("resid 0 to 184")),reference_M1)
        rmsd_M1 = np.array(rmsd_M1).reshape((len(rmsd_M1),1))
        rmsd_M2 = md.rmsd(traj.atom_slice(top.select("resid 185 to 369")),reference_M2)
        rmsd_M2 = np.array(rmsd_M2).reshape((len(rmsd_M2),1))

        ### Vector of distance between center of mass of monomer 1 and monomer 2
        centers_M1 = md.compute_center_of_mass(traj.atom_slice(top.select("resid 0 to 184")))
        centers_M2 = md.compute_center_of_mass(traj.atom_slice(top.select("resid 185 to 369")))

        V_Distance = np.array(centers_M2) - np.array(centers_M1)

        ### Three bases of monomer 1.
        ### Coordinates of three atoms chosen for deriving three bases.
        centers_M1_A = md.compute_center_of_mass(traj.atom_slice(top.select("resid 181")))
        centers_M1_B = md.compute_center_of_mass(traj.atom_slice(top.select("resid 158")))
        centers_M1_C = md.compute_center_of_mass(traj.atom_slice(top.select("resid 19")))

        ### Two vectors used for calculating bases.
        V_M1_a1 = np.array(centers_M1_B) - np.array(centers_M1_A)
        V_M1_a2 = np.array(centers_M1_C) - np.array(centers_M1_A)

        ### First base
        V_M1_x = V_M1_a1 / np.linalg.norm( V_M1_a1, axis=-1)[:, np.newaxis]

        ### Projection of a2 on V_M1_x.
        S_M1_a2_x = np.array( [np.dot(V_M1_a2[i],V_M1_x[i]) for i in range(frames)] )

        ### Second base
        V_M1_y = V_M1_a2 - np.array( [ S_M1_a2_x[i] * V_M1_x[i] for i in range(frames) ])
        V_M1_y = V_M1_y / np.linalg.norm( V_M1_y, axis=-1)[:, np.newaxis]

        ### Third base
        V_M1_z = np.cross(V_M1_x, V_M1_y)
        V_M1_z = V_M1_z / np.linalg.norm( V_M1_z, axis=-1)[:, np.newaxis]

        ### Calculate the projections of distance vector along x1, y1, z1.
        d_x     = np.array( [ [ np.dot(V_Distance[i], V_M1_x[i]) ] for i in range(frames) ] )
        d_y     = np.array( [ [ np.dot(V_Distance[i], V_M1_y[i]) ] for i in range(frames) ] )
        d_z     = np.array( [ [ np.dot(V_Distance[i], V_M1_z[i]) ] for i in range(frames) ] )

        ### Three bases of monomer 2.

        ### Coordinates of three atoms chosen for deriving three bases.
        centers_M2_A = md.compute_center_of_mass(traj.atom_slice(top.select("resid 365")))
        centers_M2_B = md.compute_center_of_mass(traj.atom_slice(top.select("resid 342")))
        centers_M2_C = md.compute_center_of_mass(traj.atom_slice(top.select("resid 203")))

        ### Two vectors used for calculating bases.
        V_M2_a1 = np.array(centers_M2_B) - np.array(centers_M2_A)
        V_M2_a2 = np.array(centers_M2_C) - np.array(centers_M2_A)

        ### First base
        V_M2_x = V_M2_a1 / np.linalg.norm( V_M2_a1, axis=-1)[:, np.newaxis]

        ### Projection of a2 on V_M2_x.
        S_M2_a2_x = np.array( [np.dot(V_M2_a2[i],V_M2_x[i]) for i in range(frames)] )

        ### Second base
        V_M2_y = V_M2_a2 - np.array( [ S_M2_a2_x[i] * V_M2_x[i] for i in range(frames) ])
        V_M2_y = V_M2_y / np.linalg.norm( V_M2_y, axis=-1)[:, np.newaxis]

        ### Third base
        V_M2_z = np.cross(V_M2_x, V_M2_y)
        V_M2_z = V_M2_z / np.linalg.norm( V_M2_z, axis=-1)[:, np.newaxis]

        ### Calculate Angles
        cos_ang_x = np.array( [ [ np.dot(V_M1_x[i], V_M2_x[i]) / (np.linalg.norm(V_M1_x[i]) * np.linalg.norm(V_M2_x[i])) ] for i in range(frames) ] )
        ang_x     = np.arccos(cos_ang_x)
        cos_ang_y = np.array( [ [ np.dot(V_M1_y[i], V_M2_y[i]) / (np.linalg.norm(V_M1_y[i]) * np.linalg.norm(V_M2_y[i])) ] for i in range(frames) ] )
        ang_y     = np.arccos(cos_ang_y)
        cos_ang_z = np.array( [ [ np.dot(V_M1_z[i], V_M2_z[i]) / (np.linalg.norm(V_M1_z[i]) * np.linalg.norm(V_M2_z[i])) ] for i in range(frames) ] )
        ang_z     = np.arccos(cos_ang_z)

        ### Save the features
        dataset   = np.hstack((d_x, d_y, d_z, ang_x, ang_y, ang_z, rmsd_M1, rmsd_M2))
        np.save( trajname + ".npy",dataset)
