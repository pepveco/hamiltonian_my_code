README – Schwinger–Dyson Integration for the p-Spin Model (Hamiltonian Dynamics)
================================================================================

This folder contains the FORTRAN codes developed for the numerical integration of the Schwinger–Dyson equations governing the Hamiltonian dynamics of the pure and mixed p-spin models.

A comprehensive description of the theoretical framework, numerical schemes, and results will be provided in the final version of Giuseppe F. Conte’s MSc thesis.

Contents
--------

The repository is organized into subfolders by model type and integration scheme. The main codes implement second-order Predictor–Corrector algorithms (with possible modifications), and are designed to simulate:

- The pure p = 2 model (integrable regime),
- The pure p = 3 model (non-integrable),
- The mixed p + s model (e.g., p = 2, s = 3).

Each implementation computes time-dependent correlation and response functions, the Lagrange multiplier enforcing the spherical constraint, and other dynamical observables (e.g., energy density, overlap functions).

Compilation
-----------

The codes are written in standard FORTRAN 90. For typical use, they can be compiled with `gfortran`:

    gfortran -g  filename.f90


For simulations involving large matrices or to optimize memory usage (especially when running on constrained RAM), we recommend compiling with flags for position-independent code and architecture-specific optimizations:

    gfortran -g -fPIC -mcmodel=medium -fPIC -o a.out filename.f90

This version is better suited for large-scale runs, as it improves performance and may reduce memory fragmentation.

Execution and Output
--------------------

The programs do not require external input files: all relevant physical parameters must be set manually inside the source code before compilation.

For the `double_corrector.f90` file (used for the pure p-spin case), the following parameters must be specified:

- `mass` — mass of the spins
- `gamma` — friction coefficient (typically zero for Hamiltonian dynamics)
- `itemp` — initila configuration temperature T'
- `J_f` — coupling constant for the post-quench dynamics (equal to `J_0` in the unquenched case)
- `p` — degree of interaction

In the case of mixed dynamics (e.g., in codes for the \( p + s \) model), the coupling constants for the initial and post-quench dynamics (`c2`, `c3`, `c2_zero`, `c3_zero`) are already hard-coded as also the degree of the interaction (p=2,s=3)

---

Upon execution, the code generates several output files containing:

- The time evolution of the Lagrange multiplier \( z(t) \), enforcing the spherical constraint;
- The two-time correlation function \( C(t, t') \) and response function \( R(t, t') \), evaluated:
  - for fixed waiting times \( t_w \) and varying \( t \);
  - or for fixed \( tw \) and varying \( t \);
- The integrated response (susceptibility) over appropriate time intervals.

Output file naming conventions may vary depending on the model and code version.
The executable runs without external input files. All parameters (e.g., time step `h`, total time, temperature `T`, initial condition settings) are defined directly in the source code.

Upon executio

--------------------

The compiled executables run without the need for input files. Simulation parameters such as time step (`h`), total simulation time, and temperature are set directly in the source code.

The output is written to plain-text files and includes:

- Correlation matrix `C(t, t')`
- Response function matrix `R(t, t')`
- Lagrange multiplier `μ(t)`
- One-time quantities and time averages

Limitations and Notes
---------------------

⚠️ **Important notes for the p = 2 case:**

- Two different subfolders are provided (`PC/` and `central_derivatoive/`) to address numerical issues specific to the integrable regime.

- At **high temperatures**, the algorithm may drift toward a non-physical asymptotic state (e.g., condensation from a paramagnetic start), due to limitations of the discretized dynamics, ONLY IN THE PC SCHEME

- At **low temperatures**, the correlation between different replicas may deviate from the true equilibrium value as a result of accumulated numerical errors during integration, IN BOTH CASES. This issue have to be solved to have a complete and unbiased understading of the results 
in a wide initial temeprature regime.

These artifacts are specific to the p = 2 case and do not affect the non-integrable p = 3 simulations.


Acknowledgments
---------------

I would like to thank Lorenzo Correale for the fruitful discussions, as well as Alessandro Tartaglia and Nicolas Nessi for providing reference codes that helped and inspired the development of this implementation, despite the fact that the approach followed here is entirely different.


License and Contact
-------------------

This code is released for academic use in connection with Giuseppe F. Conte’s MSc thesis. For questions or collaboration requests, please contact the author.

