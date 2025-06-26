– Numerical integration of Schwinger–Dyson equations for p = 2
===================================================================================

This folder contains two subdirectories, each providing a distinct implementation of the Predictor–Corrector scheme 
for the numerical integration of the Schwinger–Dyson equations in the spherical p-spin model with p = 2.

Due to numerical instabilities and structural differences in the dynamical equations for the integrable p = 2 case
(as compared to non-integrable cases such as p = 3), we developed two alternative approaches:

- `PC/` contains a direct adaptation of the Predictor–Corrector algorithm to the p = 2 setting.
- `central_derivativen/` includes an alternative version that modifies the discretization scheme and 
   the treatment of the constraint, to improve numerical stability.

These two subfolders reflect different strategies to accurately capture the dynamics in the integrable regime,
where standard methods may fail or require refinement.

⚠️ Note: Both codes are specifically designed for the p = 2 case and are not intended for use in the non-integrable case (p ≥ 3).

For more details on the algorithms, the motivation behind the two versions, and the corresponding numerical results, 
please refer to Giuseppe F. Conte’s MSc thesis.
