# PIMC1D
Path Integral MC of a single particle in 1D for training

This is a relatively simple Fortran code to simulate a single quantum particle in a 1D external potential. The Statistical Mechanis in the Canonical ensemble is in the Feynman paths representation and the statistical averages are obtained by the Metropolis Monte Carlo method of the isomorphic classical system. 

When the number of time slices is **even** the code performs a Path Integral Monte Carlo (PIMC) calculation implementing several drivers for path sampling:
1) single bead local moves (but with gaussian a priori sampling matrix) (mc);
2) Levy flight reconstruction (by successive bisections) of part of the path (levy);
3) Smart Monte Carlo moves (global) in the direct space (smc);
4) simple Metropolis MC in normal modes space with gaussian sampling matrix (fpimc);
5) Smart Monte Carlo in the normal modes space (smc_nm)

When the number of time slices is **odd** the code performs a ground state calculation with Reptation Quantum Monte Carlo (RQMC). In this case the path is open and the inverse temperature is the projection time. With a simple modification, explained in the template.in file, the end points trial function can be set to a gaussian with given variance. In the RQMC method the code implements two different drivers:
1) the strandard Reptation moves with change of growing direction at each move (reptate);
2) the bounce algorithm described in C. Pierleoni and D.M. Ceperley, ChemPhysChem 2005, 6, 1872 – 1878 (10.1002/cphc.200400587) (bounce)
