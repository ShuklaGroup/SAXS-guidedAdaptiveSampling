### Performing kinetic Monte Carlo simulations on Markov state models to estimate the sampling time required to reach the 
### target state from an arbitrary inital state using single long simulation.
### Required packages: numpy, msmbuilder
### @Chuankai Zhao, czhao37@illinois.edu

from msmbuilder.utils import io
import numpy as np

msm    =  io.load("MSM25.pkl")

### Varying the number of adaptive rounds, and the length of indidual trajectories of each round. 
Steps   =  [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
Num     =  10
Maxrun  =  100
Rounds  =  [ round + 1 for round in range(Maxrun) ]

Totraj  =  [ round*Num for round in Rounds]

### Calculate total sampling time required.
def calcSamplingtime(trajs, tot, step):
  ### lag_time in microsecond
  lag_time  = 0.05
  endframes  = [] 
  for i in range(tot):
    traj = trajs[i]
    for frame in range(step):
      state = traj[frame]
      if state in [361]:
        print("Reach folded state at traj = ", i, " frame = ", frame, "!")
        endframes.append(frame) 
        break

  if len(endframes) == 0:
    endframes.append(np.nan)
  print(endframes)

  if np.unique(endframes)[0] != np.nan:
    endframe = np.nanmin(endframes)
    endtime  = (endframe + 1) * lag_time * tot
  if np.unique(endframes)[0] == np.nan:
    endtime  = np.nan
  print(endtime)
  return endtime

### Kernel to run kinetic MC simulation
def runMCSimulation(msm,tot, step):
  trajs = []
  init  = 327
  for i in range(tot):
    traj = run(init,msm,step)
    trajs.append(traj)
  np.save("Single_MC_Traj_" + str(tot) + "_" + str(step) + ".npy", np.array(trajs))
  return trajs

def run(init,msm,step):
  try:
    traj = msm.sample_discrete(state=init, n_steps=step, random_state=None)
  except KeyError:
    print("ERROR, using random state for starting!")
    traj = msm.sample_discrete(state=None, n_steps=step)
  return traj

### Main function
SamplingTs = np.zeros(( len(Totraj), len(Steps), 3))
tot_serial = 0
for tot in Totraj:
  for step in Steps:
    print("Now, it's totrajs = ",tot, " nsteps = ", step, "!")
    trajs  = runMCSimulation(msm, tot, step)
    print(np.shape(trajs))
    time   = calcSamplingtime(trajs, tot, step)
    SamplingTs[tot_serial][step-1][0] = tot
    SamplingTs[tot_serial][step-1][1] = step
    SamplingTs[tot_serial][step-1][2] = time
  tot_serial = tot_serial + 1
np.save("SingleSampTimes.npy",SamplingTs)
