# SAXS-guided Adaptive Sampling

The pipeline of adaptive sampling consists of iteratively running short parallel simulations, clustering the trajectories based on some structural features, and seeding new simulations from certain clusters according to some selection criterion\cite{huang2009rapid}. The key of SAXS-guided adaptive sampling is to incorporate the SAXS information in the selections of seeding structures for iterative sampling. This is achieved by converting the SAXS profile into a SAXS discrepancy scoring function, which measures the degree of similarity between the target experimental or theoretical SAXS profile and the SAXS profile calculated from the structural models of each cluster. By selecting the clusters which are closer to the target, we bias the sampling direction while leaving the energy function unchanged. By iteratively running short parallel simulations, we drive the system of interest towards the target structure while still maintaining accurate thermodynamics and kinetics in the sampling process. The SAXS discrepancy function used in this study is the commonly used reduced $\chi^{2}$ function (equation \ref{eq:chi}): 
\begin{linenomath}
\begin{equation}
\chi^{2}= \frac{1}{N-1} \sum_{i=1}^{N} (\frac{\mu I_{state}(q_{i}) - I_{target}(q_{i})}{\sigma_{target} (q_{i})})^{2} 
\label{eq:chi}
\end{equation}
\end{linenomath}
where $q_{i}$ is the momentum transfer (\textit{q}=4$\pi$sin$\theta$/$\lambda$, 2$\theta$ is the scattering angle and $\lambda$ is the x-ray wave length), $I_{target}(q_{i})$ and $I_{state}(q_{i})$ are the scattering intensities of target SAXS profile and each cluster state at $q_{i}$, $\sigma_{target}(q_{i})$ is the error of target scattering intensity at $q_{i}$, \textit{N} is the total number of data points in the SAXS scattering curves, $\mu$ is a scaling factor.
