### Clustering the trajectories based on the featurized data from 01-Featurization 
### Required packages: numpy, msmbuilder
### @Chuankai Zhao, czhao37@illinois.edu

import numpy as np
import msmbuilder.cluster
from msmbuilder.utils import io

# Read the list of MD trajectories
List = [ line.rstrip() for line in open("List","r") ]

# Load the featurized data
dataset = []
for i in List:
  a = np.load( i + ".npy")
  dataset.append(a)
tran_data = np.transpose(np.vstack(dataset))

# Preprocessing: normalize the features 
dataset = []
for i in List:
  a = np.load(i + ".npy")
  b = np.array(a)
  tran_b = np.transpose(b)
  for i in range(8):
    tran_b[i] = (tran_b[i] - np.mean(tran_data[i]))/np.std(tran_data[i])
  c = np.transpose(tran_b)
  dataset.append(c)

# Perform time independent component analysis on the normalized features
from msmbuilder.decomposition import tICA
tica = tICA(n_components=4, lag_time=1)
tica.fit(dataset)
tica_traj = tica.transform(dataset)
np.save('tica_traj', tica_traj)


# Perform k-means clustering based on the first 4 tICs
states = msmbuilder.cluster.KMeans(n_clusters=200)
states.fit(tica_traj)

# Save the clustered files
io.dump(states,'clustering_tica.pkl')
