double_corrector.f90 – Predictor–Corrector integration of Schwinger–Dyson equations
====================================================================================

This folder contains the file `double_corrector.f90`, which implements a second-order 
Predictor–Corrector algorithm with two correction steps. The code numerically integrates 
the Schwinger–Dyson equations for the spherical p-spin model, with a specific focus 
on the case p = 3.

The implementation is based on a reformulated version of the dynamical equations, 
incorporating a novel treatment of the spherical constraint derived from an 
alternative analysis (the so-called "Correale trick").

⚠️ Note: This scheme is *not suitable* for the p = 2 case, as the structure of the 
equations and constraints differs substantially.

For a detailed explanation of the theoretical background, algorithmic design, 
and numerical results, please refer to Giuseppe F. Conte’s MSc thesis.

