#  Loop-extruders alter bacterial chromosome topology to direct entropic forces for segregation

This repository contains the analysis code used for the manuscript ["Loop-extruders alter bacterial chromosome topology to direct entropic forces for segregation"](https://www.biorxiv.org/content/10.1101/2023.06.30.547230v1).

There are several files in the top directory:

- Main_dynamic.jl
- Main_steady_state.jl
- fork_separation.jl
- check_time_scales.jl
- interpolate_parameters.jl
- load_specificity_and_number_smcs.jl

These files can be run to generate statistics and plots for all simulation data found in the directory "Simulation data". The large, raw simulation data used in the paper can be downloaded from Zenodo, but we provide analysed statistics in this repository.

The most time-consuming parts of the analysis are the calculations and plotting of distance and contact maps. If these maps are saved in PDF format, plotting takes many hours on a laptop. Saving the plots in png format is significantly faster, then plotting should take less than an hour. Analysis of contact and distance maps also takes hours.

Other parts of the script are much faster, and should run in the orders of minutes.

The code is written in Julia, and the required packages are listed in the file Project.toml. The easiest way to install the required packages is to use the Julia package manager. Start Julia, open Pkg mode by pressing `]`, and run the command `activate .` to activate the project environment. Then run the command `instantiate` to install the required packages.

If you have any questions, you can email j.k.harju [at] vu.nl
