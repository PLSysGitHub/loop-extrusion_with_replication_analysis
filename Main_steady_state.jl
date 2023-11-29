"""
This file runs a pipeline that:
1. Calculates statistics for saved steady state loop-extruder simulations
2. Makes and saves plots for those statistics

The script was used to generate graphs for Figure 2
"""

include("src/Steady_state_analysis.jl")
include("src/Steady_state_plots.jl")

sample_from=200 #which block number sampling should start from
repeat=false #whether or not exisiting stat files should be overwritten
monomer_size=129/1000 #μm
N=404 #number monomers
figure_type="pdf" #can also save as PDF; just change to "pdf"

#where to save data and plots
config_dir="Simulation_data/Steady_state/" #Directory where simulation output .h5 are stored
out_dir= "Stats/Steady_state_from_$sample_from/" #statistics are saved here
plot_dir= "Plots/Steady_state_from_$sample_from/" #plots are saved here
comparison_plot_dir= "Plots/Steady_state_compare_from_$sample_from/" #plots comparing simulations with and without SMCs
convergence_dir="Stats/Steady_state_convergence/" #stats as a function of simulation time
convergence_plot_dir="Plots/Steady_state_convergence/" #plots to check for convergence

#first do mean z-positions, segregated fractions
calc_fast_steady_state_stats(config_dir,out_dir, skip_done=!repeat, sample_from=sample_from)

#Hi-C and distance maps are slow to calculate; comment if not needed 
calc_slow_steady_state_stats(config_dir, out_dir, skip_done=!repeat, sample_from=sample_from)
make_slow_eq_plots(out_dir, plot_dir, monomer_size, filetype=figure_type)

#Compare steady state data from simulations with or without SMCs; make plots

no_smcs=["Steady_state_No_smcs_no_tether/","Steady_state_No_smcs/", "Steady_state_No_smcs_no_tether/"]
ideal=["Steady_state_Ideal_no_tether/","Steady_state_Ideal/", "Steady_state_Ideal_no_tether_tied_forks/"]
smcs=["Steady_state_Nontopological_smcs_no_tether/", "Steady_state_Nontopological_smcs/", "Steady_state_No_smcs_no_tether_tied_forks/"]
lbls=["No_tether/", "Tethers/", "Tied_forks/"]
Rs=[[100,200, 250, 300,402], [100, 200, 300, 402], [100, 200, 300, 402]]
for (index,no_smcs_dir) in enumerate(no_smcs)
    smcs_dir=smcs[index]
    ideal_dir=ideal[index]
    lbl=lbls[index]
    if index==3
        compare_point_data(Rs[index], out_dir*no_smcs_dir, out_dir*smcs_dir, out_dir*ideal_dir,
            comparison_plot_dir*lbl, monomer_size, N, filetype=figure_type, lbl2="Tied forks", lbl1="Free forks", lbl3="Tied forks, Ideal", line2=:dashdot)
    else
        compare_point_data(Rs[index], out_dir*no_smcs_dir, out_dir*smcs_dir, out_dir*ideal_dir,
            comparison_plot_dir*lbl, monomer_size, N, filetype=figure_type)
    end
end

#Make plots that track stats as a function of simulation time; check for convergence

check_convergence(config_dir,convergence_dir)
make_convergence_plots(convergence_dir,convergence_plot_dir, monomer_size, filetype=figure_type)