p2_xentral_derivative_no_corr.f90 – Numerical integration of Schwinger–Dyson equations for p = 2 using central derivative discretization
=================================================================================================================================================
The file `p2_central_derivative_no_corr.f90` implements a central-difference integration scheme for the case p = 2,
 designed to overcome the numerical instabilities encountered (expecially) at high temperatures when using the Predictor–Corrector (PC) algorithm.
 However, we have numerical deviations from the right behavior in the low temperature regime. 
 These issues are specific to this setting and not related to the formulation proposed by Lorenzo Correale for p = 3.

