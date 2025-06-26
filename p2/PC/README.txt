p2_double_corrector.f90 – Numerical integration of Schwinger–Dyson equations for p = 2 using PC algorithm
=============================================================================================================

This folder contains two subdirectories, each providing a distinct implementation of the Predictor–Corrector (PC) 
scheme for the numerical integration of the Schwinger–Dyson equations in the spherical p-spin model with p = 2.

Due to differences in the structure of the equations compared to the p = 3 case, the standard integration scheme encounters several issues:

- At **high temperatures**, the code may converge toward an incorrect asymptotic state. 
  Specifically, starting from a paramagnetic initial condition, the system may drift toward a condensed state
   that is not consistent with the true equilibrium solution.
  
- At low temperatures, the algorithm exhibits oscillations that progressively deviate the value of the correlation 
between different replicas from the true asymptotic equilibrium value, due to the accumulation of numerical errors during the integration.


For a detailed explanation of the theoretical motivation, the behavior across temperature regimes, 
and the derivation of the additional equations, please refer to Giuseppe F. Conte’s MSc thesis.
