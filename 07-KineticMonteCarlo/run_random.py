### Performing kinetic Monte Carlo simulations on Markov state models to estimate the sampling time required to reach the
### target state from an arbitrary inital state using random adaptive sampling strategy.
### Required packages: numpy, msmbuilder
### @Chuankai Zhao, czhao37@illinois.edu

from msmbuilder.utils import io
import numpy as np

msm    =  io.load("MSM25.pkl")

Steps   =  [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
Num     =  10
Maxrun  =  100
Rounds  =  [ round + 1 for round in range(Maxrun) ]

def runMCSimulation(msm, round, step, Num):
  init   = 327                        # unfolded state
  trajs  = []

  R=round
  initRoundFlag = True
  if initRoundFlag:
    for i in range(Num):
      traj      = run(init, msm, step)
      traj      = [ traj ]
      trajs.append(traj)
    initRoundFlag = False

  unique = np.unique(trajs)
  inits  = np.zeros(Num)
  for i in range(Num):
    inits[i]  = int(np.random.choice(unique,1))
  
  if round > 1 and not initRoundFlag:
    round = round - 1
    for iter in range(round):
      inits_iter = inits 
      for i in range(Num):
        traj     = run(inits_iter[i], msm, step)
        trajs[i].append(traj) 

      unique   = np.unique(trajs)
      inits    = np.zeros(Num)
      for i in range(Num):
        inits[i]  = int(np.random.choice(unique,1))

  np.save("random_MC_Traj_" + str(R) + "_" + str(step) + "_" + str(Num) + ".npy", np.array(trajs))
  return trajs

def calcSamplingtime(trajs, round, step, Num):
  lag_time  = 0.05
  endframes  = [] 
  for i in range(Num):
    notreachflag = True
    if notreachflag:
      for iter in range(round):
        traj = trajs[i][iter]
        for frame in range(step):
          state = traj[frame]
          if state in [361]:
            tot_frame = iter * step + frame + 1
            endframes.append(tot_frame)
            print("Reach folded states at i = ", i, " iter = ", iter, " frame = ", frame)
            notreachflag = False
            break
        if not notreachflag:
          break

  if len(endframes) == 0:
    endframes.append(np.nan)

  print(endframes)

  if np.unique(endframes)[0] != np.nan:
    endframe = np.nanmin(endframes)
    endtime  = endframe * lag_time * Num
  if np.unique(endframes)[0] == np.nan:
    endtime  = np.nan

  print(endtime)
  return endtime

def run(init,msm,step):
  try:
    traj = msm.sample_discrete(state=init, n_steps=step, random_state=None)
  except KeyError:
    print("ERROR, using random state for starting!")
    traj = msm.sample_discrete(state=None, n_steps=step)
  return traj

SamplingTs = np.zeros(( len(Rounds), len(Steps), 3))
for round in Rounds:
  for step in Steps:
    print("Now, it's round = ",round, " nsteps = ", step, "!") 
    trajs  = runMCSimulation(msm, round, step, Num)
    time   = calcSamplingtime(trajs, round, step, Num)
    SamplingTs[round-1][step-1][0] = round
    SamplingTs[round-1][step-1][1] = step
    SamplingTs[round-1][step-1][2] = time
np.save("RandomSampTimes.npy",SamplingTs)
