# PANTHER
Physics-based semi-ANalytical Tool for Human-induced Earthquake Rupture

Model framework for stress changes and fault reactivation due to reservoir pressure and temperature changes

- 2D plane-strain model
- two reservoir comparments offset by a fault
- pressure and/or temperature changes in one or both compartments
- pressure and/or temperature diffusion to the seal and base
- poro-elastic and thermo-elastic stress changes based on Jansen et al.
- fault reactivation and aseismic slip
- stochastic analysis
- parallel computing

How to run?
- clone the repo git@github.com:TNO/PANTHER.git
- open the Panther.prj file in Matlab
- run panther.m or check out the various example files

Requires:
- Matlab
- Parallel Processing Toolbox (optional)
