### Given the SAXS profiles of the target and individual states, calculate the SAXS discrepancy scores
### @Chuankai Zhao, czhao37@illinois.edu

import numpy as np
import math
import glob

### load the target SAXS profiles, including both intensities and errors
f_saxs_data = np.loadtxt("avg_native.dat")
f_saxs = np.transpose(f_saxs_data)[1]
f_saxs_err = np.transpose(f_saxs_data)[2]

N = 51
scores = []

### calculate the reduced chi^2 SAXS discrepancy scores between each state and the target 
for i in range(500):
  state_score = []
  for j in range(100):
    s_saxs = np.loadtxt("STATE" + str(i) + "_" + str(j) + ".txt")
    s_saxs = np.transpose(s_saxs)[1]
    #s_saxs = (f_saxs[0]/s_saxs[0])*s_saxs
    sum    = 0.
    for k in range(N):
      sum = sum + ( (f_saxs[k] - s_saxs[k])/f_saxs_err[k] )**2
    sum    = sum/(N-1)
    state_score.append(sum)
  scores.append([i, np.mean(state_score), np.std(state_score)])

### save the SAXS discrepancy scores 
np.savetxt("ProteinG_discrepancy.txt", np.array(scores))
