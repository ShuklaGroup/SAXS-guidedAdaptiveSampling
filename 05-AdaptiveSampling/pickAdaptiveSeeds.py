### Given the SAXS discrepancy scores, pick the N states with the least SAXS discrepancy scores
### for next round of adaptive sampling. Other information can be incorporated together with SAXS
### to pick the seeding structures for adaptive sampling. 
### Required packages: numpy
### Output: input files for CPPTRAJ to extract the seeding structures and save as rst (Amber restart files) 

import numpy as np

### read the list of MD trajectories
List = [ line.rstrip() for line in open("List", "r") ]

### read the saxs discrepancy scores and sort based on scores
saxs = np.loadtxt("ProteinG_discrepancy.txt")
saxs.view('i8,i8,i8').sort(order=['f1'], axis=0)

### load the cluster files
cl = io.load('../cluster/clustering_tica.pkl')
cl = cl.labels_

### randomly pick M structures from each of the N states with the lowest SAXS discrepancy scores for adaptive sampling.
M=5
N=10

import random

### define the list of frames to extract from the MD trajectories
final_selects = []

for i in range(N):
  selects = []
  for j in range(len(cl)):
    for k in range(len(cl[j])):
      if cl[i][j] == saxs[i][0]:
        selects.append[j][k]
  sel = np.random.choice(range(len(selects)),size=M)
  for j in range(M):
    final_selects.append(sel[j])


### write the CPPTRAJ input file to extract frames
round   = 1
path    = "/home/czhao37/3-ABA/trajs_full/"
opath   = "/home/czhao37/3-ABA/Round" + str(round+1) + "/0-Minimization/"
selects = final_selects

for i in range(len(selects)):
    select = selects[i]
    select = [ int(select[0]), int(select[1]) ]
    print(select)
    name = "r" + str(round+1) + "par" + str(i+1)
    f = open("cpptraj_" + name + ".in", "w")
    f.write("parm " + path + List[select[0]] + ".prmtop\n")
    f.write("trajin " + path + List[select[0]] + ".mdcrd " + str(select[1]+1) + " " + str(select[1]+2) +  " 4\n")
    f.write("autoimage\n")
    f.write("parmbox alpha 90 beta 90 gamma 90\n")
    f.write("trajout " + opath + name + ".rst restart\n")
    f.write("parmwrite out " + opath + name + ".prmtop\n")
