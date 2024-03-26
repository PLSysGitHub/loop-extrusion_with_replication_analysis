"""
This file contains scripts for plots comparing simulation results with different rates of specific loading and with different numbers of loop-extruders.
"""

include("src/Plot_functions.jl")
using LaTeXStrings

no_smcs_file_name(R)= "Stats/Steady_state_from_200/Steady_state_No_smcs_no_tether/R_$(R)_N_404_colRate_0.03_stoch_0.05_sps_2500_trunc_1.5_pullF_2.0_top_mon_0/steady_state_stats.h5"
smcs_file_name(R, parS, stallFork)= "Stats/Steady_state_from_200/Steady_state_Nontopological_smcs_no_tether/R_$(R)_stallFork_$(stallFork)_N_404_colRate_0.03_sep_101_parSstrength_$(parS)_lifetime_40400_stoch_0.05_sps_2500_trunc_1.5_pullF_2.0/steady_state_stats.h5"
dynamic_smcs_file_name(parS, stallFork)= "Stats/Segregation_dynamic/Nontopological_smcs_no_tether/stallFork_$(stallFork)_N_404_colRate_0.03_sep_101_parSstrength_$(parS)_lifetime_40400_rateFork_0.07440000000000001_stoch_0.05_sps_2500_trunc_1.5_pullF_2.0_segregation.h5"
dynamic_no_smcs_file_name= "Stats/Segregation_dynamic/No_smcs_no_tether/N_404_colRate_0.03_rateFork_0.07440000000000001_stoch_0.05_sps_2500_trunc_1.5_pullF_2.0_segregation.h5"

"""
Given a set of parS strengths, fetch the steady state segregated fraction and its standard deviation for each parS strength.
"""
function fetch_steady_state_seg_indices_vary_parS_strength(parS_strengths, R=402, stallFork=0.0)
    seg_fractions=[]
    seg_fraction_stds=[]
    
    for parS_strength in parS_strengths
        file_name=smcs_file_name(R, parS_strength, stallFork)
        h5open(file_name, "r") do file
            seg_fraction=read(file, "seg_index")
            seg_fraction_std=read(file, "seg_index_std")
    
            push!(seg_fractions, seg_fraction)
            push!(seg_fraction_stds, seg_fraction_std)
        end
    end

    return seg_fractions, seg_fraction_stds
end

"""
Given a set of parS strengths, fetch the steady state segregated fraction and its standard deviation for each parS strength.
"""
function fetch_steady_state_mean_zs_seg_frac_for_parS_strength(parS_strength, R=402, stallFork=1.0)    
    file_name=smcs_file_name(R, parS_strength, stallFork)
    return h5read(file_name, "mean_zs"), h5read(file_name, "seg_index"), h5read(file_name, "seg_index_std")
end

"""
Given a set of parS strengths, fetch the segregated fraction as a function of replicated length and its standard deviation for each parS strength.
"""
function fetch_dynamic_seg_indices_vary_parS_strength(parS_strengths, stallFork=0.0)

    seg_fraction_trajectories=[]
    seg_fraction_stds_trajectories=[]
    
    for parS_strength in parS_strengths
        file_name=dynamic_smcs_file_name(parS_strength, stallFork)
    
        h5open(file_name, "r") do file
            seg_fraction=read(file, "seg_index")
            push!(seg_fraction_trajectories, seg_fraction)
            seg_fraction_std=read(file, "seg_index_std")
            push!(seg_fraction_stds_trajectories, seg_fraction_std)
        end
    end

    return seg_fraction_trajectories, seg_fraction_stds_trajectories
end

"""
Given a set of separations, fetch the dynamic segregated fraction and the number of loop-extruders for each separation.
"""
function fetch_dynamic_seg_indices_vary_SMC_number(separations, N=4040)
    numbers=[]
    seg_fractions=[]
    
    for sep in separations
        push!(numbers, floor(Int,N/sep))
        file_name="Stats/Segregation_dynamic/Nontopological_smcs_no_tether/stallFork_0.0_N_404_colRate_0.03_sep_$(sep)_parSstrength_4040.0_lifetime_40400_rateFork_0.07440000000000001_stoch_0.05_sps_2500_trunc_1.5_pullF_2.0_segregation.h5"
        h5open(file_name, "r") do file
            seg_fraction=read(file, "seg_index")
            push!(seg_fractions, seg_fraction)
        end
    end

    return seg_fractions, numbers
end


#Plotting directory
plot_dir="Plots/Compare_offloading_and_num_smcs/"
if !isdir(plot_dir)
    mkpath(plot_dir)
end

#For x-axis ticks
N=404
r_Mb=collect(0:N)./100

#First, steady state with varying levels of loading specificity
parS_strengths=["1010.0", "2020.0", "4040.0", "40400.0", "404000.0"] #parS strengths to compare

seg_fractions, seg_fraction_stds=fetch_steady_state_seg_indices_vary_parS_strength(parS_strengths)

no_smcs_file=no_smcs_file_name(402)
no_smcs_segregated_fraction=h5read(no_smcs_file, "seg_index")
no_smcs_segregated_fraction_std=h5read(no_smcs_file, "seg_index_std")

#Plot of segregated fraction at R=N as a function of parS strength
plot(parse.(Float64, parS_strengths), seg_fractions, yerr=seg_fraction_stds, xlabel=L"\mathrm{Relative\ loading\ affinity\ } ori", 
    ylabel=L"\mathrm{Mean\ segregated\ fraction,\ }R=402", color=:black, ylims=(0,1), label="Loop-extruders", xscale=:log10, size=(700,400))
hline!([no_smcs_segregated_fraction], ribbon=no_smcs_segregated_fraction_std, color=:teal, label="No loop-extruders")
savefig(plot_dir*"Steady_state_segregation_vs_parS_strength.pdf")

seg_fracs=[]
seg_frac_std=[]
seg_fracs_no_smcs=[]
seg_frac_std_no_smcs=[]
#Compare the segregated fractions at different Rs for highest affinity
Rs=[100, 200, 300, 402]
for R in Rs
    monomer_size=129/1000
    mean_zs, seg_f, s_std=fetch_steady_state_mean_zs_seg_frac_for_parS_strength(parS_strengths[end], R)
    push!(seg_fracs, seg_f)
    push!(seg_frac_std, s_std)
    push!(seg_fracs_no_smcs, h5read(no_smcs_file_name(R), "seg_index"))
    push!(seg_frac_std_no_smcs, h5read(no_smcs_file_name(R), "seg_index_std"))
end

plot_compare_segregated_fractions(Rs, hcat(seg_fracs_no_smcs, seg_fracs),hcat(seg_frac_std_no_smcs,seg_frac_std))
savefig(plot_dir*"compare_segregated_fractions_vs_R_affinity_$(parS_strengths[end]).pdf")

#Now the dynamic case. Compare segregated fraction over time.
parS_strengths=[4040.0, 40400.0, 404000.0] #parS strengths to compare
labels=["Affinity ori=$(round(Int,af))" for af in parS_strengths]

seg_fraction_trajectories, seg_fraction_stds_trajectories=fetch_dynamic_seg_indices_vary_parS_strength(parS_strengths)
seg_fraction_trajectory_no_smcs=h5read(dynamic_no_smcs_file_name, "seg_index")

alphas=[1-i*(0.9)/length(parS_strengths) for i in 0:length(parS_strengths)-1]
plot(r_Mb, seg_fraction_trajectories, color=:black, alpha=reshape(alphas,1,:), ylims=(0,1),
    xlabel="Replicated length [Mb]", ylabel="Segregated fraction", label=reshape(labels,1,:), size=(800,400))
plot!(r_Mb, seg_fraction_trajectory_no_smcs, color=:teal, label="No loop-extruders")
savefig(plot_dir*"Dynamic_segregation_trajectories.pdf")

#Compare with or without off-loading at replication forks, targeted loading
seg_offload=h5read(dynamic_smcs_file_name(404000.0, 1.0), "seg_index")
seg_dont_offload=h5read(dynamic_smcs_file_name(404000.0, 0.0), "seg_index")

plot(r_Mb, seg_offload, color=:black, label=L"P_U=1", ylabel="Segregated fraction", xlabel="Replicated length [Mb]", size=(800,400), ylims=(0,1))
plot!(r_Mb, seg_dont_offload, color=:gray, label=L"P_U=0")
plot!(r_Mb, seg_fraction_trajectory_no_smcs, color=:teal, label="No loop-extruders")
savefig(plot_dir*"off-loading_changes_affinity_404000.pdf")

#Comparing different numbers of loop-extruders
separations=[50, 101, 202, 808]

seg_fractions, numbers=fetch_dynamic_seg_indices_vary_SMC_number(separations)
colors=[:black for i in 1:length(numbers)]
alphas=[1-(0.95/length(numbers))*i for i in 0:length(numbers)-1]

#add the no loop-extruders case
push!(seg_fractions, seg_fraction_trajectory_no_smcs)
push!(numbers, 0)
push!(colors, :teal) #same color as in main figures
push!(alphas, 1.0)

#Plot all on same axes
plot(r_Mb, seg_fractions, label=reshape(numbers, 1,:), xlabel="Replicated length [Mb]", ylabel="Segregated fraction", size=(650,400), ylims=(0,1),color=reshape(colors, 1,:), alpha=reshape(alphas,1,:))
savefig(plot_dir*"seg_index_vary_num_smcs.pdf")