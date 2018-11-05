# SAXS-guided Adaptive Sampling

Molecular dynamics (MD) simulations can be utilized to predict protein structure ensembles and dynamics, though sufficient sampling of molecular ensembles and identification of key biologically relevant conformations remains challenging. SAXS-guided adaptive sampling is an efficient unbiased sampling protocol for large-scale MD simulations to enhance time efficiency and identify functionally relevant conformations of proteins and complexes. The central idea is to utilize protein structural information obtained from low-resolution experimental techniques to guide Markov state model (MSM)-based adaptive sampling. 

## MSM-based adaptive sampling

The pipline of MSM-based adaptive sampling consists of (1) iteratively running short parallel simulations, (2) clustering the trajectories based on some structural features, and (3) seeding new simulations from certain clusters according to some selection criterion. For example, iteratively picking the clusters with the least population is one of the common strategies in sampling protein dynamics. 


## Incorporation of SAXS information in MSM-based adaptive sampling
Small angle x-ray scattering (SAXS) is popular technique to characterize both structured and intrinsically disordered biomolecules in solution, especially for complexes. The key of SAXS-guided adaptive sampling is to incorporate the SAXS information in the selections of seeding structures for iterative sampling. This is achieved by converting the SAXS profile into a SAXS discrepancy scoring function, which measures the degree of similarity between the target experimental or theoretical SAXS profile and the SAXS profile calculated from the structural models of each cluster. By selecting the clusters which are closer to the target, we bias the sampling direction while leaving the energy function unchanged. By iteratively running short parallel simulations, we drive the system of interest towards the target structure while still maintaining accurate thermodynamics and kinetics in the sampling process. 

## Implementation of SAXS-guided adaptive sampling
This collection includes the python scripts for implementing SAXS-guided adaptive sampling. More details for each individual step are described below. 

### 01-Featurization
Beginning with running short parallel simulations on available protein structures (e.g. crystal structures), we obtain a collection of protein conformations. In order to perform clustering on these conformations, we need to calculate some structural features that can describe protein conformational changes and dynamics, known as the featurization process. 

In this folder, three sample scripts demonstrate the featurizations based on root mean square deviation (RMSD) of proteins from native structure (featurizeRMSD.py) and dihedral angles (featurizeDihedral.py) for protein folding trajetories, as well as a combination of various features (featurizeAssociation.py) for protein-protein association process. 

### 02-Clustering
Unsupervised clustering can then be performed on the obtained features using a variety of algorithms. Time-lagged independent component (tICA) analysis is commonly used in clustering MD trajectories, in order to capture the slowest-relaxing degrees of freedom from linear combinations of input features. In this folder, a sample script (clustering_tica.py) demonstrates the clustering process using tICA and k-means algorithm with MSMBuilder. 

### 03-ExtractStructures
Next, in order to calculate the SAXS profiles for each cluster, we need to extract protein structures from each cluster as PDB or trajectory format depending on the choice of software for SAXS calculation. 

In this folder, the first set of scripts (0-ExtractToPDB) show the processes of randomly extracting structures from each cluster as PDB files, and using the PDB files to calculate the SAXS profiles using Crysol. Crysol is the most commonly used implicit-solvent SAXS modeling software. Another set of scripts (1-ExtractToXTC) show the processes of extracting all frames from each cluster in MD data (with water molecules) as single xtc trajectories. The xtc tracjetories can then be used as input files to calculate the SAXS profiles using WAXSiS. WAXSiS is an explicit-solvent SAXS modeling software. 

### 04-SAXSDiscrepancy
The script (discrepancy.py) is used to calculate the SAXS discrepancy values between the SAXS profiles of individual states and the target experimental or computation SAXS data. The SAXS discrepancy score implemented is the reduced chi^2 function. 

### 05-AdaptiveSampling
The calculated SAXS discrepancy scores are then used to choose the clusters for starting new MD simulations. The script (pickAdaptiveSeeds.py) is used to randomly pick a few structures from the clusters (with the lowest SAXS discrepancy scores) as the adaptive seeds for next round of parallel simulations. 

Additional structural information such as distance restraints inferred from evolutionary couplings can be combined with SAXS to pick the structures, in order to collectively enhance sampling efficiency. Steps 01-05 are repeated until the target structure is discovered. 

### 06-MSMConstructions
For the analysis, in order to construct the MSM to describe protein dynamics, one will have to choose the lag time. The script msmTimescales.py is used to plot the timescale plots with respect to lag time. By checking the convergence of the slowest timescales, one can justify the choice of lag time.  

### 07-KineticMonteCarlo
In this folder, the set of python scripts are used to perform kinetic Monte Carlo simulations on the constructed MSMs to compare the sampling efficiency using a variety of adaptive sampling protocols. 

## Reference
Please refer to the following manuscript for more method details.

Zhao, C. & Shukla, D. SAXS-guided enhanced unbiased sampling for structure determination of proteins and complexes. Sci. Rep. (2018). Under review.  
