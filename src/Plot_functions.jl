"""
Given that statistics have been calculated, this file can be used to make the plots
"""

using Plots, DelimitedFiles, HDF5
include("helpers.jl")

pythonplot(grid=false,label="",framestyle=:box,colorbar=true, linewidth=2,
    guidefontsize=15, tickfontsize=15,colorbar_tickfontsize=12,legend=:outerright, markersize=7,
    colorbar_titlefontsize=15, legendfontsize=15, xticks=(1:100:404, string.([2, 3, 0, 1, 2])),
    palette=:tab10)

"""
Either pass a Hi-C maps for a single chromosome, or a replicated one (set replicated=true)

Plots the Hi-C map with monomer mid at the center of the plot,
if R is not NaN, draw dashed lines to indicate fork positions
"""
function plot_hic(M; R=NaN, mid=1, replicated=false)
    N=size(M)[1]
    cticks=0:0.02:0.1

    if replicated
        M_shifted=shifted_map_replicated(M,mid)
    else
        M_shifted=shifted_map(M,mid)
    end
    pl=heatmap(M_shifted,clims=(0,0.1),color=cgrad(:BuPu), size=(500,450),
            xlims=(1,N),ylims=(1,N), xlabel="Genomic position [Mb]",
            colorbar_ticks=(cticks,cticks), aspect_ratio=1., 
            ylabel="Genomic position [Mb]", colorbartitle="Contact frequency",
            yticks=(100:100:N,[3, 0, 1, 2]))
    if !isnan(R)
        pl=vline!(pl,[N/2+R/2, N/2-R/2], color=:red, linestyle=:dash)
        pl=hline!(pl,[N/2+R/2, N/2-R/2], color=:red, linestyle=:dash)
    end
end

"""
Given separated Hi-C maps for the two strands, and the trans-contact map,
plots the Hi-C map with monomer mid at the center of the plot,
if R is not NaN, draw dashed lines to indicate fork positions
"""
function plot_hic(old, new, inter; R=NaN, mid=1)
        N=size(old)[1]
        cticks=0:0.02:0.1

        all_hic=zeros(2*N, 2*N)
    
        all_hic[1:N,1:N].=old
        all_hic[N+1:2*N, N+1:2*N].=new
        all_hic[1:N, N+1:2*N].=transpose(inter)
        all_hic[N+1:2*N, 1:N].=inter

        M_shifted=shifted_map_replicated(all_hic,mid)

        pl=heatmap(M_shifted,clims=(0,0.1),color=cgrad(:BuPu), size=(500,450),
                xlims=(1,2*N),ylims=(1,2*N),
                colorbar_ticks=(cticks,cticks), aspect_ratio=1.,
                ylabel="Strand 1 / 2 [Mb]", colorbartitle="Contact frequency",
                xticks=(vcat(1:100:301, 501:100:801), ["2", "3", "0", "1", "3'", "0'", "1'", "2'"]), 
                yticks=(vcat(100:100:301, 501:100:801), ["3", "0", "1", "3'", "0'", "1'", "2'"]),
                xlabel="Strand 1 / 2 [Mb]")
        if !isnan(R)
            pl=vline!(pl,[N/2+R/2, N/2-R/2, 3*N/2+R/2, 3*N/2-R/2], color=:red, linestyle=:dash)
            pl=hline!(pl,[N/2+R/2, N/2-R/2, 3*N/2+R/2, 3*N/2-R/2], color=:red, linestyle=:dash)
        end
end

"""
Given inter-chromosomal contacts, plots the trans-contact map with monomer mid at the center of the plot,
if R is not NaN, draw dashed lines to indicate fork positions
"""
function plot_inter_hic(M; R=NaN, mid=1)
        N=size(M)[1]
        cticks=0:0.02:0.1

        M_shifted=shifted_map(M,mid)

        pl=heatmap(M_shifted,clims=(0,0.1),color=cgrad(:BuPu), size=(500,450),
                xlims=(1,N),ylims=(1,N), xlabel="Genomic position, Strand 1 [Mb]",
                colorbar_ticks=(cticks,cticks), aspect_ratio=1.,
                ylabel="Genomic position, Strand 2 [Mb]", colorbartitle="Contact frequency",
                yticks=(100:100:N,["3'", "0'", "1'", "2'"]))
        if !isnan(R)
                pl=vline!(pl,[N/2+R/2, N/2-R/2], color=:red, linestyle=:dash)
                pl=hline!(pl,[N/2+R/2, N/2-R/2], color=:red, linestyle=:dash)
        end
end

########################
# Average distances

"""
Given a full distance map (i.e. with both strands), plots the average distance map with monomer mid at the center of the plot,
draw dashed lines to indicate fork positions

Set scale_μm to the monomer size in μm
"""
function plot_distance_matrix(M, fork; scale_μm=1)
    d_lim=3.2
    N=Int(size(M)[1]/2)
    M_copy=copy(M)*scale_μm

    #unreplicated regions to max value
    M_copy[fork[1]+N:fork[2]+N,:].=d_lim
    M_copy[:,fork[1]+N:fork[2]+N].=d_lim
    
    #shift the origin to mid-strand
    M_copy=shifted_map_replicated(M_copy)

    heatmap(M_copy,color=cgrad(:thermal, rev=true),clims=(0,d_lim),
            aspect_ratio=1,xlims=(1,2*N),ylims=(1,2*N),
            xticks=(vcat(1:100:301, 501:100:801), ["2", "3", "0", "1", "3'", "0'", "1'", "2'"]), 
            yticks=(vcat(100:100:301, 501:100:801), ["3", "0", "1", "3'", "0'", "1'", "2'"]),
            xlabel="Strand 1 / 2 [Mb]", ylabel="Strand 1 / 2 [Mb]", colorbartitle="Average distance [μm]")
        vline!(vcat(fork.+[N/2,-N/2], fork.+[3*N/2,N/2]), color=:red, linestyle=:dash)
        hline!(vcat(fork.+[N/2,-N/2], fork.+[3*N/2,N/2]), color=:red, linestyle=:dash)
end

"""
Plot the inter-chromosomal average distance map
set scale_μm to the monomer size in μm
draw dashed lines to indicate fork positions
"""
function plot_distance_matrix_oldnew(M, fork, scale_μm=1)
    d_lim=3.2
    N=Int(size(M)[1]/2)
    M_copy=M[N+1:end,1:N]*scale_μm
    M_copy[fork[1]:fork[2],:].=d_lim

    heatmap(shifted_map(M_copy,1),color=cgrad(:thermal, rev=true),clims=(0,d_lim),
            aspect_ratio=1,xlims=(1,N),ylims=(1,N),
            yticks=(100:100:N,["3'", "0'", "1'", "2'"]),
            xlabel="Genomic position, Strand 1 [Mb]", ylabel="Genomic position, Strand 2 [Mb]", colorbartitle="Average distance [μm]")

    vline!([N/2+fork[1], fork[2]-N/2], color=:red, linestyle=:dash)
    hline!([N/2+fork[1], fork[2]-N/2], color=:red, linestyle=:dash)
end


"""
Given average long axis positions from simulations with and without smcs,
plot them on the same figure, with an indication of the cell size and fork positions
"""
function plot_compare_av_z(av1, av2, fork, N=404)
        ter=ceil(Int,N/2)

        unrepl1=fork[1]:ter-1
        unrepl2=ter:fork[2]

        repl=vcat(fork[2]:N, 1:fork[1])
        repl_inds= repl.-(ter)
        repl_inds[repl_inds.<=0].+=N+1

        L=fork_to_height(fork)

        old=av2[1:N]
        new=av2[N+1:end]

        hline([L/2], ribbon=(L, 0), color=:turquoise, fillalpha=0.2)
        hline!([-L/2], color=:turquoise)
        hline!([0], color=:grey, linealpha=0.3)
        plot!(repl_inds,old[repl], ylabel="Long axis position [μm]",xlabel="Genomic position [Mb]", color=1)
        plot!(unrepl2.-(ter),old[unrepl2], color=:black,label="Loop-extruders")
        plot!(unrepl1.-(ter-N+1),old[unrepl1], color=:black)
        plot!(repl_inds,new[repl], color=2, size=[700,400])

        old2=av1[1:N]
        new2=av1[N+1:end]

        plot!(repl_inds,old2[repl], color=1, linestyle=:dash, ylims=(-3.45/2, 3.45/2))
        plot!(unrepl2.-(ter),old2[unrepl2], color=:black,label="No loop-extruders", linestyle=:dash)
        plot!(unrepl1.-(ter-N),old2[unrepl1], color=:black,linestyle=:dash)
        plot!(repl_inds,new2[repl], color=2, linestyle=:dash)
end

"""
Given average long axis positions from simulations at different time points,
plot them on the same figure
"""
function plot_av_zs_all(avs, forks, N=405)
        ter=ceil(Int,N/2)
        nt=size(avs)[2]

        hline([0], color=:grey, linealpha=0.3, ylims=(-1.5, 1.5))

        for i in 1:nt
                α=1-(i-1)/(nt)

                unrepl1=forks[1,i]:ter-1
                unrepl2=ter:min(N,forks[2,i])

                repl=vcat(min(N,forks[2,i]-1):N, 1:forks[1,i]+1)
                repl_inds= repl.-(ter)
                repl_inds[repl_inds.<=0].+=N+1

                old=avs[1:N,i]
                new=avs[N+1:end,i]
                t=Int(round(fork_to_time(forks[:,i],N=N), digits=-1))
                if i==1
                        unrepl2=vcat(ter:N,1:ter)
                        plot!(old[unrepl2], color=:black, label="0 min")
                else 
                        plot!(repl_inds,old[repl], ylabel="Long axis position [μm]",xlabel="Genomic position [Mb]", color=1, alpha=α)
                        plot!(repl_inds,new[repl], color=2, size=[550,400], alpha=α)
                        plot!(unrepl2.-(ter),old[unrepl2], color=:black, alpha=α, label="$t min")
                        plot!(unrepl1.-(ter-N+1),old[unrepl1], color=:black, alpha=α)
                end

        end
end


"""
Plot a kymograph of the origins and terminus positions along the long axis as a function of time
Include an indicator for the cell size as a function of time
"""
function ori_ter_kymograph(zs, stds)
        N=size(zs)[2]-1
        hs=map(i->R_to_height(i), 0:N)

        ter=floor(Int,N/2)

        plot(hs./2, ribbon=(hs, zeros(N+1)), color=:turquoise, fillalpha=0.2)
        plot!(-hs./2, color=:turquoise)
        plot!(zs[1,:], ribbon=stds[1,:], label="Ori 1", color=1, ylabel="Mean long axis position")
        plot!(zs[N+1,:], ribbon=stds[N+1,:], label="Ori 2", color=2, xlabel="Replicated length [Mb]")
        plot!(zs[ter,:], ribbon=stds[ter,:], label="Ter", color=:black, alpha=0.8, xticks=(1:100:404, string.([0,1, 2, 3, 4])))
end

"""
Plot a kymograph of the long axis positions of the centers of mass of the chromosomal strands as a function of time
Include an indicator for the cell size as a function of time
"""
function com_kymograph(zs, stds)
        N=size(zs)[2]-1
        hs=map(i->R_to_height(i), 0:N)

        plot(hs./2, ribbon=(hs, zeros(N+1)), color=:turquoise, fillalpha=0.2)
        plot!(zs[1,:], ribbon=stds[1,:], label="Strand 1", color=1, ylabel="Center of mass\nMean long axis position")
        plot!(zs[2,:], ribbon=stds[2,:], label="Strand 2", color=2, xlabel="Replicated length [Mb]")
        plot!(zs[3,:], ribbon=stds[3,:], label="Ter strand", color=:black, alpha=0.8, xticks=(1:100:404, string.([0,1, 2, 3, 4])))
        plot!(-hs./2, color=:turquoise)
end

"""
Given a set of replicated lengths and segregated fractions, plot a scatter of the fractions
with error bars

segregated_fractions should be an array of size (2, n) where n is the number of replicated lengths
The first row should be the segregated fractions without loop-extruders, the second row with loop-extruders
"""
function plot_compare_segregated_fractions(Rs, segregated_fractions, segregated_fractions_std)
        scatter(Rs, segregated_fractions[2,:], ribbon=segregated_fractions_std[2,:], fillalpha=0.2,
                ylabel="Segregated fraction", xlabel="Replicated length [Mb]", marker=:circle, markersize=10,size=(800,400),
                xticks=(100:100:404, string.([1, 2, 3, 4])), label="Loop-extruders", ylims=(0,1), color=:black)
        scatter!(Rs, segregated_fractions[1,:], ribbon=segregated_fractions_std[1,:], fillalpha=0.4,
                label="No loop-extruders", color=:grey, marker=:diamond, markersize=10)
end

"""
Given long vectors of segregated fractions and their standard deviations as a function of R from dynamic simulations
Plot the segregated fractions as a function of R

segregated_fractions_1 should correspond to data without loop-extruders, segregated_fractions_2 with loop-extruders
"""
function plot_compare_segregated_fractions(segregated_fractions_1, segregated_fractions_std_1, segregated_fractions_2, segregated_fractions_std_2)
        plot(10:N, segregated_fractions_2, ribbon=segregated_fractions_std_2, ylims=(0,1), xlabel="Replicated length [Mb]", ylabel="Segregated fraction",
                xticks=(1:100:404, string.([0, 1, 2, 3, 4])), label="Loop-extruders", color=:black, size=(800,400), fillalpha=0.2)
        plot!(10:N, segregated_fractions_1, ribbon=segregated_fractions_std_1, label="No loop-extruders", linestyle=:dash, color=:grey, fillalpha=0.4)
end

"""
Given long vectors of the radial separation between the centers of mass of replicated strands, plot the radial separation as a function of R

radial_distances_1 should correspond to data without loop-extruders, radial_distances_2 with loop-extruders
"""
function plot_compare_radial_distances(radial_distances_1, radial_distances_std_1, radial_distances_2, radial_distances_std_2)
        plot(radial_distances_2, ribbon=radial_distances_std_2, label="Loop-extruders", ylabel="Radial separation\ncenters of mass [μm]", size=(800,400),
                xticks=(1:100:404, string.([0, 1, 2, 3, 4])),xlabel="Replicated length [Mb]", color=:black, fillalpha=0.2)
        plot!(radial_distances_1, ribbon=radial_distances_std_1, label="No loop-extruders", linestyle=:dash, color=:grey, fillalpha=0.4)
end

"""
Given long vectors of the long axis positions of monomers as a function of simulation time, plot the long axis positions
        of the oris and the ter
"""
function plot_convergence_oris(long_axis_positions, R)
        L=R_to_height(R)
        L_max=R_to_height(404)*1.1

        hline([L/2], ribbon=(L, 0), color=:turquoise, fillalpha=0.2)
        hline!([-L/2], color=:turquoise)
        plot!(long_axis_positions[1,:], xlabel="Simulation time", ylabel="Average long axis position [μm]", label="Ori 1", color=1)
        plot!(long_axis_positions[N+1,:], label="Ori 2", ylims=(-L_max/2,L_max/2), legend=:outerright, color=2)
        plot!(long_axis_positions[Int(N/2),:], label="Ter", xticks=:auto, color=:black)
end