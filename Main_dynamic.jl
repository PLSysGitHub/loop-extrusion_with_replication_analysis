"""
This file runs a pipeline that:
1. Calculates statistics for saved replicating chromosome simulations
2. Makes and saves plots for those statistics
3. Creates an example animation for each simulation, as well as individual snapshot .pngs

The script was used to generate plots for Figures 4, S2, and S3
"""

include("src/Dynamic_analysis.jl")
include("src/Dynamic_plots.jl")

repeat=true #should existing statistics files be overwritten?
times=round.(Int,R_to_time.([0, 1000, 2020, 3000, 4020])) #Which times (min) do we (roughly) want samples for?
N=404 #total number monomers
monomer_size=129/1000 #μm

# Directory names
data_parent_dir="Simulation_data/Dynamic/" #Directory where simulation output .h5 are stored
out_parent_dir="Stats/Dynamic/"
plot_dir="Plots/Dynamic/"
comparison_plot_dir="Plots/Dynamic_comparison/"
out_segr_dir="Stats/Segregation_dynamic/"
plots_segr_dir="Plots/Segregation_dynamic/"
animation_dir="Animations/Dynamic/"
snapshot_dir="Snapshots/Dynamic/"


#################################
#          Run pipeline         #
#################################

#Set skip_done parameter to false to overwrite old results.
calc_time_point_stats(times, N, data_parent_dir,out_parent_dir, skip_done=!repeat)
make_time_point_plots(times, out_parent_dir, plot_dir,monomer_size,N=N)

#Segregation plots
calc_time_course_stats(N, data_parent_dir,out_segr_dir, skip_done=!repeat)
make_time_course_plots(out_segr_dir,plots_segr_dir, monomer_size, N)

# Plots to compare with or without smcs

#list the files you want to compare to each other
no_smcs=["No_smcs/N_404_colRate_0.03_rateFork_0.07440000000000001_stoch_0.05_sps_2500_trunc_1.5_pullF_2.0",
    "No_smcs_no_tether/N_404_colRate_0.03_rateFork_0.07440000000000001_stoch_0.05_sps_2500_trunc_1.5_pullF_2.0",
    "No_smcs/N_404_colRate_0.03_rateFork_0.07440000000000001_stoch_0.05_sps_2500_trunc_1.5_pullF_2.0",
    "No_smcs_no_tether/N_404_colRate_0.03_rateFork_0.07440000000000001_stoch_0.05_sps_2500_trunc_1.5_pullF_2.0",
    "No_smcs_no_tether/N_404_colRate_0.03_rateFork_0.07440000000000001_stoch_0.05_sps_2500_trunc_3.0_pullF_2.0",
    "No_smcs_no_tether/N_404_colRate_0.03_rateFork_0.07440000000000001_stoch_0.05_sps_250_trunc_1.5_pullF_2.0",
    "No_smcs_no_tether/N_404_colRate_0.3_rateFork_0.07440000000000001_stoch_0.05_sps_2500_trunc_1.5_pullF_2.0"
    ]

smcs=["Nontopological_smcs/stallFork_0.0_N_404_colRate_0.03_sep_101_parSstrength_4040.0_lifetime_40400_rateFork_0.07440000000000001_stoch_0.05_sps_2500_trunc_1.5_pullF_2.0",
    "Nontopological_smcs_no_tether/stallFork_0.0_N_404_colRate_0.03_sep_101_parSstrength_4040.0_lifetime_40400_rateFork_0.07440000000000001_stoch_0.05_sps_2500_trunc_1.5_pullF_2.0",
    "Topological_smcs_no_tether/stallFork_0.0_N_404_colRate_0.03_sep_101_parSstrength_4040.0_lifetime_40400_rateFork_0.07440000000000001_stoch_0.05_sps_2500_trunc_1.5_pullF_2.0",
    "Topological_smcs_no_tether/stallFork_0.0_N_404_colRate_0.03_sep_101_parSstrength_4040.0_lifetime_40400_rateFork_0.07440000000000001_stoch_0.05_sps_2500_trunc_1.5_pullF_2.0",
    "Nontopological_smcs_no_tether/stallFork_0.0_N_404_colRate_0.03_sep_101_parSstrength_4040.0_lifetime_40400_rateFork_0.07440000000000001_stoch_0.05_sps_2500_trunc_3.0_pullF_2.0",
    "Nontopological_smcs_no_tether/stallFork_0.0_N_404_colRate_0.03_sep_101_parSstrength_4040.0_lifetime_40400_rateFork_0.07440000000000001_stoch_0.05_sps_250_trunc_1.5_pullF_2.0",
    "Nontopological_smcs_no_tether/stallFork_0.0_N_404_colRate_0.3_sep_101_parSstrength_4040.0_lifetime_40400_rateFork_0.07440000000000001_stoch_0.05_sps_2500_trunc_1.5_pullF_2.0"
    ]

labels=["Nontopological_tether/", "Nontopological_no_tether/","Topological_tether/",  "Topological_no_tether/", "Nontopological_no_tether_3_trunc/", "Nontopological_no_tether_sps_250/", "Non_topological_no_tether_colRate_0.3/"] #names for subdirectories where comparisons saved

for (index, no_smc_file) in enumerate(no_smcs)
    smc_file=smcs[index]
    lbl=labels[index]
    
    plot_compare_time_point_data(times,out_parent_dir*no_smc_file*".h5", out_parent_dir*smc_file*".h5", comparison_plot_dir*lbl, monomer_size,N)
    plot_compare_segregation_data(out_segr_dir*no_smc_file*"_segregation.h5", out_segr_dir*smc_file*"_segregation.h5",comparison_plot_dir*lbl,monomer_size,N)
end

#Save some animations of trajectories
include("src/Trajectory_animations.jl")
for which in 1:50:200 #which simulation file to animate
    one_animation_per_sim(data_parent_dir,animation_dir; skip_done=!repeat, which=which)
end
snapshots_per_sim(data_parent_dir, snapshot_dir; which=101)