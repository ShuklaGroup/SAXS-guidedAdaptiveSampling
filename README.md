# SAXS-guided Adaptive Sampling

Molecular simulations can be utilized to predict protein structure ensembles and dynamics, though sufficient sampling of molecular ensembles and identification of key biologically relevant conformations remains challenging. SAXS-guided adaptive sampling is an efficient sampling protocol for large-scale unbiased molecular dynamics (MD) simulations to enhance time efficiency and identify functionally relevant conformations of proteins and complexes. The central idea of this methodology is to utilize coarse protein structural information obtained from low-resolution experimental techniques to guide Markov state model (MSM)-based adaptive sampling. 

## MSM-based adaptive sampling

The pipline of MSM-based adaptive sampling consists of 1) iteratively running short parallel simulations, 2) clustering the trajectories based on some structural features, and 3) seeding new simulations from certain clusters according to some selection criterion. 


## Incorporation of SAXS information in adaptive sampling
The key of SAXS-guided adaptive sampling is to incorporate the SAXS information in the selections of seeding structures for iterative sampling. This is achieved by converting the SAXS profile into a SAXS discrepancy scoring function, which measures the degree of similarity between the target experimental or theoretical SAXS profile and the SAXS profile calculated from the structural models of each cluster. By selecting the clusters which are closer to the target, we bias the sampling direction while leaving the energy function unchanged. By iteratively running short parallel simulations, we drive the system of interest towards the target structure while still maintaining accurate thermodynamics and kinetics in the sampling process. 

