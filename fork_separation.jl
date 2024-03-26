"""
This script was used to compare fork separations from steady state simulations
to expectation values found by assuming a linear dependence of strand extension 
on strand length.
"""

include("src/helpers.jl")
using HDF5, Plots, StatsBase

pythonplot(grid=false,label="",framestyle=:box,colorbar=true, linewidth=1.5,
    guidefontsize=15, tickfontsize=15,colorbar_tickfontsize=12,legend=:outerright, markersize=7,
    colorbar_titlefontsize=15, legendfontsize=15)

"""
Assume that all chromosome segments fill roughly same fraction of long axis

Then segment between forks, of length n=min(R,N-R), fills L * n/(N/2),
where L is available long axis distance, and N is total number of monomers

"""
function predicted_fork_sep(R,N,L)
    if R>=N/2
        return L*(N-R)/R #Circle of length 2R spread across length L, linear segment N-R
    else
        return L*R/(N/2) #Circle of length N spread across length L, linear segment R
    end
end

"""
Assume that all chromosome segments fill roughly same fraction of long axis

With smcs, ori-ori-ter organization then implies that the distance from the ori to forks is proportional to R/n
After R>N/2, ori-ter-ori organization means forks at mid-cell
"""
function predicted_fork_height_smcs(R,N,L)
    if R>=N/2
        z=0 #middle of the linearized largest circle
    else
        z= L*R/N-L/2 #circle of length N, genomic length R/2 from ori to fork. Mid-cell at zero.
    end
    return z
end

"""
Go through the steady state results directory, extract the mean fork positions and separations for all found values of R
"""
function fetch_fork_positions_separations_no_smcs(stats_directory)
    #Collect the mean fork separations and standard deviations from steady state sims
    Rs=[0] #replicated lengths read from file names
    seps=[0.] #fork separations
    stds=[0.] #separation standard deviation

    fpos=[[0.,0.]] #fork positions
    fstds=[[0.,0.]] #position standard deviation

    #Read simulation data files for long axis separation of forks
    for dir in subdirs(stats_directory*"/Steady_state_No_smcs_no_tether/")
        R=parse_after(dir)
        push!(Rs, R)
        h5open(dir*"/steady_state_stats.h5", "r") do f
            separations=read(f, "z_dist")
            std=read(f, "z_dist_std")

            push!(seps, separations[end])
            push!(stds, std[end])

            z=read(f, "mean_zs")
            std=read(f, "zs_std")
            fork=read(f, "fork")
            push!(fpos, z[fork])
            push!(fstds, std[fork])
        end
    end

    seps=seps[sortperm(Rs)]
    stds=stds[sortperm(Rs)]
    fpos=fpos[sortperm(Rs)]
    fstds=fstds[sortperm(Rs)]
    sort!(Rs)

    #convert fpos and fstds to arrays
    fpos=reduce(hcat, fpos)
    fstds=reduce(hcat, fstds)

    fpos[:,1].=NaN
    return Rs, seps, stds, fpos, fstds
end

"""
Go through the steady state results directory, extract the mean fork positions and separations for all found values of R
"""
function fetch_fork_positions_separations_with_smcs(stats_directory)
    Rs=[0]
    zs=[[predicted_fork_height_smcs(0,N,R_to_height(0)-diameter)/monomer_size,predicted_fork_height_smcs(0,N,R_to_height(0)-diameter)/monomer_size]]
    zstds=[[0.,0.]]

    fseps=[0.]
    fstds=[0.]

    #Read simulation data files for long axis separation of forks
    for dir in subdirs(stats_directory*"Steady_state_Nontopological_smcs_no_tether/")
        if contains(dir, "lifetime_40400_") && contains(dir, "parSstrength_4040.0_")
            R=parse_after(dir)
            push!(Rs, R)

            h5open(dir*"/steady_state_stats.h5", "r") do f
                z=read(f, "mean_zs")
                std=read(f, "zs_std")
                fork=read(f, "fork")
                push!(zs, z[fork])
                push!(zstds, std[fork])

                separations=read(f, "z_dist")
                std=read(f, "z_dist_std")

                push!(fseps, separations[end])
                push!(fstds, std[end])

            end
        end
    end

    zs=zs[sortperm(Rs)]
    zstds=zstds[sortperm(Rs)]
    fseps=fseps[sortperm(Rs)]
    fstds=fstds[sortperm(Rs)]
    sort!(Rs)

    zs=reduce(hcat, zs)
    zstds=reduce(hcat, zstds)

    return Rs, zs, zstds, fseps, fstds
end

#Directory names
stats_directory="Stats/Steady_state_from_200/"
if !isdir(stats_directory)
    @error "Could not find statistics for comparison! Please run Main_steady_state.jl first!"
end

plot_directory="Plots/Steady_state_compare_from_200/"
if !isdir(plot_directory)
    mkpath(plot_directory)
end

fig_type=".pdf" #can also save png or other formats

#Parameters for cell dimensions
diameter=6.4*129/1000 #μm
monomer_size=129/1000 #μm
N=404 #number of monomers

#fetch data for no smcs
Rs, seps, stds, fpos, fstds=fetch_fork_positions_separations_no_smcs(stats_directory)

Ls=R_to_height.(0:N)
predicted_fork_separations=predicted_fork_sep.(0:N,N,Ls.-diameter)

#Now plot the mean fork positions
plot(95:N, Ls[95:N]./2, ribbon=(Ls[95:N], zeros(N+1)), color=:turquoise, fillalpha=0.2, size=(530,400), xticks=(100:100:400, 1:4));
plot!(95:N, -Ls[95:N]/2, color=:turquoise, xlims=(95,N), xlabel="Replicated length [Mb]", ylabel="Mean fork long axis pos. [μm]");
scatter!(Rs, fpos'*monomer_size, ribbon=fstds'*monomer_size, markersize=5, color=[:red :lightcoral], label=["Fork 1" "Fork 2"])
savefig(plot_directory*"fork_position_no_LE"*fig_type)

#Fork heights in the presence of loop-extruders can be predicted similarly
Rs, zs, zstds, fseps, fstds=fetch_fork_positions_separations_with_smcs(stats_directory)

#Calculate predictions, assuming that available space is L-d (fork in middle of confinement blob)
predicted_fork_pos=predicted_fork_height_smcs.(0:N,N,Ls.-diameter)

#Mean long-axis positions of forks with SMCs
plot(0:N, Ls./2, ribbon=(Ls, zeros(N+1)), color=:turquoise, fillalpha=0.2, size=(530,400));
plot!(0:N, -Ls/2, color=:turquoise, xlims=(95,N));
scatter!(Rs, zs'*monomer_size, ribbon=zstds'*monomer_size, markersize=5, color=[:red :lightcoral], label=["Fork 1" "Fork 2"], ylabel="Mean fork long axis pos. [μm]");
plot!(0:N,predicted_fork_pos, color=:black, xlabel="Replicated length [Mb]", xticks=(100:100:400, 1:4))
savefig(plot_directory*"fork_position_LE"*fig_type)

#Fork separations, with or without loop-extruders
plot(0:N, Ls, ribbon=(Ls, zeros(N+1)), color=:turquoise, fillalpha=0.2, size=(600,400), ylabel="Mean fork long axis dist. [μm]");
scatter!(Rs, fseps*monomer_size, ribbon=fstds*monomer_size, markersize=5, color=:black, label="Loop-extruders")
plot!(0:N,predicted_fork_separations, color=:black, xlabel="Replicated length [Mb]", xticks=(0:100:400, 0:4))
scatter!(Rs, seps*monomer_size, ribbon=stds*monomer_size, markersize=5, color=:teal, label="No loop-extruders")
savefig(plot_directory*"fork_separations"*fig_type)