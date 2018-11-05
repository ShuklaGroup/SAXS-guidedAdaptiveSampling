### After running splitTrajs.py, using ./merge.sh to concanenate the block trajectories with the same cluster index
### Output: xtc trajectories including all the frames from each state.
### Output is used as the input to calculate the SAXS profiles of individual clustersi (e.g. using WAXSiS).
### Usage: ./merge.sh
### @Chuankai Zhao, czhao37@illinois.edu

#!/bin/bash

### Path to raw trajectories
path=/home/czhao37/2-SAXS-Adaptive_Samping/0-ABA-Dimerization/Round03/trajs

### Rearrange splitted small trajectories into subdirectories based on the cluster index
for i in `seq 0 199`
do
  mkdir STATE$i
  mv *STATE_${i}_* STATE$i/
  echo "STATE${i}" >> stateList
done

mkdir mergedtrajs

### Concanenate the block trajectories with the same cluster index
while read line
do
  cd ${path}/${line}
  ls -l *STATE_* > temp.txt
  awk '{ print $9 }' temp.txt > trajList
  rm -rf temp.txt
  trajNUM=`cat trajList | wc -l`
  echo $trajNUM
  if [ "$trajNUM" -eq "1" ]
  then
      mv *STATE_*.xtc ${line}.xtc
  else
      cp ../ABA_Dimer_SetUp.pdb .
      python ../trajmerge.py -i trajList -t ABA_Dimer_SetUp.pdb -o ${line}.xtc
  fi
  mv ${line}.xtc ../mergedtrajs/
done < stateList
