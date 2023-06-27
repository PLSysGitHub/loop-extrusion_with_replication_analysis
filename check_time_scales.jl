"""
This file was used to check whether the diffusion of the chromosome segments in the simulations
occurred at the expected rate. The script was used to generate Fig. FIXME

The diffusion constant for the origin in the presence of loop-extruders was of the same order of magnitude
as that found by Weber et.al. 2012. However, the average over all monomers was an order of magnitude higher.

"""

using HDF5, Plots, StatsBase, LinearAlgebra, CurveFit,LaTeXStrings,LinearAlgebra
include("src/helpers.jl")
pythonplot(grid=false,label="",framestyle=:box,colorbar=true, linewidth=1.5,
    guidefontsize=15, tickfontsize=15,colorbar_tickfontsize=12,legend=:outerright, markersize=7,
    colorbar_titlefontsize=15, legendfontsize=15)

#Constants
monomer_size=129/1000 #μm
nt=203
N=404

ori_data=[]
labels=[]

for subdir in subdirs("Simulation_data/Initial_stages/Initial_Nontopological_smcs_no_tether/")
    initial_stage_files=readdir(subdir, join=true)
    plot_directory=replace(subdir,"Simulation_data/Initial_stages/"=>"Plots/Diffusion/")*"/"
    if !isdir(plot_directory)
        mkpath(plot_directory)
        
        #Make arrays and set parameters
        num_files=length(initial_stage_files)
        long_axis_r_squared=zeros(N,nt)
        short_axis_r_squared=zeros(N,nt)
        r_squared=zeros(N,nt)

        num_samples=zeros(nt)

        for (index, f) in enumerate(initial_stage_files)
            h5open(f) do file
                ter_pos=zeros(3,N,nt)
                lowest_block=100000
                highest_block=0
                try
                    for i in keys(file)
                        block=read_attribute(file[i], "block")+1
                        h_factor=time_to_height(0)/time_to_height(block_to_time(block))

                        lowest_block=min(lowest_block, block)
                        highest_block=max(highest_block, block)

                        ter=fetch_pos_at_ind(file, i, 2, no_turn=true, N=N, spacer=1, ter_distinct=true)[:,1:N]
                        #ter[3]*=h_factor
                        ter_pos[:,:,block].=ter.*monomer_size
                    end
                    for (index,i) in enumerate(lowest_block:highest_block)
                        for j in 1:N
                            r_squared[j,index]+=sum((ter_pos[:,j,i].-ter_pos[:,j,lowest_block]).^2)
                            long_axis_r_squared[j,index]+=sum((ter_pos[3,j,i]-ter_pos[3,j,lowest_block]).^2)
                            short_axis_r_squared[j,index]+=sum((ter_pos[1,j,i]-ter_pos[1,j,lowest_block]).^2)
                        end
                        num_samples[index]+=1
                    end
                catch
                    @warn "Could not read file $f"
                end
            end
        end

        #Divide by number of samples
        for i in 1:N
            r_squared[i,:]./=num_samples
            long_axis_r_squared[i,:]./=num_samples
            short_axis_r_squared[i,:]./=num_samples
        end

        #Fit power law to mean squared distance for each monomer
        vals=map(i->power_fit(block_to_time(2:30)*6, r_squared[i,2:30]), 1:N)
        Ds=map(i->vals[i][1], 1:N)./4
        γs=map(i->vals[i][2], 1:N)

        histogram(log10.(Ds), xlabel="Diffusion coefficient D", ylabel="Number of monomers", xticks=([-2.5, -2, -1.5], [L"10^{-2.5}", L"10^{-2}", L"10^{-1.5}"]), color=:grey)
        png(plot_directory*"diffusion_early_D.png")

        histogram(γs, xlabel="Exponent, γ", ylabel="Number of monomers",color=:grey, xticks=:auto)
        png(plot_directory*"diffusion_early_γ.png")

        #Look at mean mean squared distance over all monomers and fit power law
        av_R2=transpose(mean(r_squared, dims=1))
        av_D, av_γ=power_fit(block_to_time(2:30)*6, av_R2[2:30])
        av_D/=4

        #Plot mean squared distance over all monomers
        plot(block_to_time(2:100)*6, av_R2[2:100], scale=:log10, xlabel="Simulation time, [s]", 
            ylabel="MSD monomer average, [μm²]", ylims=(10^-3,1), color=:black)
        plot!(1:0.1:1000,t-> 4*av_D*t^av_γ, color=:grey, label="Fit", alpha=0.4, xticks=:auto)
        png(plot_directory*"diffusion_early_mean_over_monomers.png")

        #Plot mean squared distance for origin
        plot(block_to_time(2:100)*6, r_squared[1,2:100], scale=:log10, xlabel="Simulation time, [s]", 
            ylabel="MSD ori, [μm²]", ylims=(10^-3,1), xlims=(1,1200), color=:black)
        plot!(1:0.1:1000,t-> 4*Ds[1]*t^γs[1], color=:grey, alpha=0.4, label="Fit", xticks=:auto)
        png(plot_directory*"diffusion_early_ori.png")

        push!(ori_data, r_squared[1,2:100])
        push!(labels, replace(subdir,"Simulation_data/Initial_stages/Initial_Nontopological_smcs_no_tether/"=>""))

        println("For $subdir, got ori parameters D=$(Ds[1]), γ=$(γs[1])")
    end
end

labels[contains.(labels,"trunc_3.0_")].="Excl. vol. x 2"
labels[contains.(labels,"colRate_0.3_")].= "Drag x 10"
labels[contains.(labels, "sps_250_")].= "3D steps / 10"
labels[contains.(labels, "stallFork")].= "Used parameters"

ori_data_array=(hcat(ori_data...))

plot(block_to_time(2:100)*6, ori_data_array, scale=:log10, xlabel="Simulation time, [s]", ticks=:auto, size=(700,400),
    ylabel="MSD ori, [μm²]", xlims=(1,1200), ylims=(10^-3,1), palette=:RdBu_4, label=reshape(labels, (1,length(labels))))
png("Plots/Diffusion/diffusion_early_ori_compare")