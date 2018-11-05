### Performing kinetic Monte Carlo simulations on Markov state models to estimate the sampling time required to reach the
### target state from an arbitrary inital state using adaptive sampling guided by SAXS and evolutionary couplings.
### Required packages: numpy, msmbuilder
### @Chuankai Zhao, czhao37@illinois.edu

from msmbuilder.utils import io
import numpy as np

msm    =  io.load("MSM25.pkl")

Steps   =  [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
Num     =  10
Maxrun  =  100
Rounds  =  [ round + 1 for round in range(Maxrun) ]

def loadSAXSdata():
  file      = open("ProteinG_discrepancy.txt")
  SAXS      = np.loadtxt(file)
  N_states  = np.shape(SAXS)[0]
  SAXS_dict = {}
  for i in range(N_states):
    state     = int(SAXS[i][0])
    intensity = SAXS[i][1]
    SAXS_dict[state] = intensity
  return SAXS_dict

def loadECdata():
  EC      = np.load("ProteinG_EC.npy")
  N_states  = np.shape(EC)[0]
  EC_dict = {}
  for i in range(N_states):
    state     = int(EC[i][0])
    dist      = EC[i][1]
    EC_dict[state] = dist
  return EC_dict

def runSX_EC_MCSimulation(msm, round, step, Num, SAXS_dict, EC_dict):
  init   = 327                        # unfoled state
  trajs  = []

  R=round
  initRoundFlag = True
  inits     = np.zeros(Num)
  trajs_rnd = [] 
  if initRoundFlag:
    for i in range(Num):
      traj      = run(init, msm, step)
      trajs_rnd.append(traj)
      traj      = [ traj ]
      trajs.append(traj)
    unique    = np.unique(trajs_rnd)
    SAXS_diff = [ SAXS_dict.get(state) for state in unique ]

    if np.shape(unique)[0] < 10:
      for i in range(np.shape(unique)[0]):
        inits[i] = int(unique[i])
      for i in range(10-np.shape(unique)[0]):
        inits[np.shape(unique)[0] + i] = int(np.random.choice(unique,1))
 
    if np.shape(unique)[0] == 10:
      for i in range(np.shape(unique)[0]):
        inits[i] = int(unique[i])
    
    if 20 > np.shape(unique)[0] > 10:
      sorted_SAXS_diff = np.unique(SAXS_diff)
      for i in range(Num):
        index = 0
        for saxs in SAXS_diff:
          if saxs == sorted_SAXS_diff[i]:
            inits[i] = int(unique[index])
          index = index + 1

    if np.shape(unique)[0] >= 20:
      selects = []
      sorted_SAXS_diff = np.unique(SAXS_diff)
      for i in range(int(2*Num)):
        index = 0
        for saxs in SAXS_diff:
          if saxs == sorted_SAXS_diff[i]:
            selects.append(int(unique[index]))
          index = index + 1
      EC_diff   = [ EC_dict.get(state) for state in selects ]
      sorted_EC_diff = np.unique(EC_diff)
      for i in range(Num):
        index = 0
        for ec in EC_diff:
          if ec == sorted_EC_diff[i]:
            inits[i] = int(selects[index])
          index = index + 1

    initRoundFlag = False
  
  if round > 1 and not initRoundFlag:
    round = round - 1
    for iter in range(round):
      inits_iter = inits 
      inits      = np.zeros(Num)
      trajs_rnd  = []
      for i in range(Num):
        traj     = run(inits_iter[i], msm, step)
        trajs[i].append(traj)
        trajs_rnd.append(traj)

      unique    = np.unique(trajs_rnd)
      SAXS_diff = [ SAXS_dict.get(state) for state in unique ]

      if np.shape(unique)[0] < 10:
        for i in range(np.shape(unique)[0]):
          inits[i] = int(unique[i])
        for i in range(10-np.shape(unique)[0]):
          inits[np.shape(unique)[0] + i] = int(np.random.choice(unique,1))
 
      if np.shape(unique)[0] == 10:
        for i in range(np.shape(unique)[0]):
          inits[i] = int(unique[i])

      if 20 > np.shape(unique)[0] > 10:
        sorted_SAXS_diff = np.unique(SAXS_diff)
        for i in range(Num):
          index = 0
          for saxs in SAXS_diff:
            if saxs == sorted_SAXS_diff[i]:
              inits[i] = int(unique[index])
            index = index + 1

      if np.shape(unique)[0] >= 20:
        selects = []
        sorted_SAXS_diff = np.unique(SAXS_diff)
        for i in range(int(2*Num)):
          index = 0
          for saxs in SAXS_diff:
            if saxs == sorted_SAXS_diff[i]:
              selects.append(int(unique[index]))
            index = index + 1
        EC_diff   = [ EC_dict.get(state) for state in selects ]
        sorted_EC_diff = np.unique(EC_diff)
        for i in range(Num):
          index = 0
          for ec in EC_diff:
            if ec == sorted_EC_diff[i]:
              inits[i] = int(selects[index])
            index = index + 1

  np.save("SAXS_EC_MC_Traj_" + str(R) + "_" + str(step) + "_" + str(Num) + ".npy", np.array(trajs))
  return trajs

def calcSamplingtime(trajs, round, step, Num):
  lag_time  = 0.05
  endframes  = [] 
  for i in range(Num):
    for iter in range(round):
      traj = trajs[i][iter]
      for frame in range(step):
        state = traj[frame]
        if state in [361]:
          tot_frame = iter * step + frame + 1
          endframes.append(tot_frame)
          break
      if len(endframes) >= i + 1: break

  if len(endframes) == 0:
    endframes.append(-np.inf)

  print(endframes)

  endframe = np.min(endframes)
  endtime  = endframe * lag_time * Num
  return endtime

def run(init,msm,step):
  try:
    traj = msm.sample_discrete(state=init, n_steps=step, random_state=None)
  except KeyError:
    print("ERROR, using random state for starting!")
    traj = msm.sample_discrete(state=None, n_steps=step)
  return traj

### Main Program ###

SAXS_dict  = loadSAXSdata()
EC_dict    = loadECdata()
SamplingTs = np.zeros(( len(Rounds), len(Steps), 3))
for round in Rounds:
  for step in Steps:
    trajs  = runSX_EC_MCSimulation(msm, round, step, Num, SAXS_dict, EC_dict)
    print(np.shape(trajs))
    time   = calcSamplingtime(trajs, round, step, Num)
    SamplingTs[round-1][step-1][0] = round
    SamplingTs[round-1][step-1][1] = step
    SamplingTs[round-1][step-1][2] = time
np.save("SAXS_EC_SamplingTs.npy",SamplingTs)
