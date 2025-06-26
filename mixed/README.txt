mixed_double_corrector.f90 – Predictor–Corrector integration of Schwinger–Dyson equations for Hamiltonina dynamics of the mixed spin model
=========================================================================================================================================

This folder contains the file `mixed_double_corrector.f90`, which implements a second-order 
Predictor–Corrector algorithm with two correction steps. The code numerically integrates 
the Schwinger–Dyson equations for the spherical mixed (p+s)-spin model, with a specific focus 
on the case p = 2 and s=3. 

The implementation is based on a reformulated version of the dynamical equations, 
incorporating a novel treatment of the spherical constraint derived from an 
alternative analysis (the so-called "Correale trick").

In order to correctly probe the integrable dynamics of the p = 2 model in the presence 
of a non-integrable term arising from the s = 3 contribution, we weight the respective 
terms associated with each part of the dynamics by appropriate coefficients. 
These coefficients tune the relative importance of each interaction term, 
effectively incorporating the coexistence of the mixed dynamics.

⚠️ Note: This scheme is *not suitable* for the p = 2 case, as the structure of the 
equations and constraints differs substantially.

For a detailed explanation of the theoretical background, algorithmic design, 
and numerical results, please refer to Giuseppe F. Conte’s MSc thesis.

