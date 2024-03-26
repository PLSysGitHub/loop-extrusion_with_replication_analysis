"""
Given that statistics have been calculated, this file can be used to make the plots
"""

using Plots, DelimitedFiles, HDF5
include("helpers.jl")

pythonplot(grid=false,label="",framestyle=:box,colorbar=true, linewidth=2,
    guidefontsize=15, tickfontsize=15,colorbar_tickfontsize=12,legend=:outerright, markersize=7,
    colorbar_titlefontsize=15, legendfontsize=15,palette=:tab10)


"""
Helper function that returns xticks in Mb
"""
function xticks_Mb(N, ori_center=true, total_Mb=4.04)
        num_Mb=floor(Int, total_Mb)
        ter=floor(Int, total_Mb/2)
        increment=round(Int, N/num_Mb)
        tick_pos=0:increment:N

        if ori_center
                tick_labels=string.(vcat(ter:num_Mb-1, 0:ter))
        else
                tick_labels=string.(0:num_Mb)
        end
        return (tick_pos, tick_labels)
end


"""
Helper function that returns xticks in Mb
"""
function yticks_Mb(N, ori_center=true, total_Mb=4.04)
        num_Mb=floor(Int, total_Mb)
        ter=floor(Int, total_Mb/2)
        increment=round(Int, N/num_Mb)
        tick_pos=increment:increment:N

        if ori_center
                tick_labels=string.(vcat(ter+1:num_Mb-1, 0:ter))
        else
                tick_labels=string.(1:num_Mb)
        end
        return (tick_pos, tick_labels)
end


"""
Helper function that returns xticks in Mb, for two chromosomes
"""
function xticks_Mb_two_chromosomes(N, ori_center=true, total_Mb=4.04)
        num_Mb=floor(Int, total_Mb)
        ter=floor(Int, num_Mb/2)
        increment=round(Int, N/num_Mb)
        tick_pos=0:increment:N
        tick_pos=vcat(tick_pos[1:end-1], tick_pos[2:end].+N)

        if ori_center
                tick_labels=string.(vcat(ter:num_Mb-1, 0:ter-1))
                tick_labels_2=string.(vcat(ter+1:num_Mb-1, 0:ter)).*"'"
                tick_labels=vcat(tick_labels, tick_labels_2)
        else
                tick_labels=string.(0:num_Mb-1)
                tick_labels_2=string.(1:num_Mb).*"'"
                tick_labels=vcat(tick_labels, tick_labels_2)
        end
        return (tick_pos, tick_labels)
end

"""
Helper function that returns xticks in Mb, for two chromosomes
"""
function yticks_Mb_two_chromosomes(N, ori_center=true, total_Mb=4.04)
        num_Mb=floor(Int, total_Mb)
        ter=floor(Int, num_Mb/2)
        increment=round(Int, N/num_Mb)
        tick_pos=0:increment:N
        tick_pos=vcat(tick_pos[2:end-1], tick_pos[2:end].+N)

        if ori_center
                tick_labels=string.(vcat(ter+1:num_Mb-1, 0:ter-1))
                tick_labels_2=string.(vcat(ter+1:num_Mb-1, 0:ter)).*"'"
                tick_labels=vcat(tick_labels, tick_labels_2)
        else
                tick_labels=string.(1:num_Mb-1)
                tick_labels_2=string.(1:num_Mb).*"'"
                tick_labels=vcat(tick_labels, tick_labels_2)
        end
        return (tick_pos, tick_labels)
end


"""
Either pass a Hi-C maps for a single chromosome, or a replicated one (set replicated=true)

Plots the Hi-C map with monomer mid at the center of the plot,
if R is not NaN, draw dashed lines to indicate fork positions
"""
function plot_hic(M; R=NaN, mid=1, replicated=false, colorlimit=0.1)
    N=size(M)[1]
    cticks=0:0.02:0.1

    if replicated
        M_shifted=shifted_map_replicated(M,mid)
    else
        M_shifted=shifted_map(M,mid)
    end
    pl=heatmap(M_shifted,clims=(0,colorlimit),color=cgrad(:BuPu), size=(500,450),
            xlims=(1,N),ylims=(1,N), xlabel="Genomic position [Mb]",
            colorbar_ticks=(cticks,cticks), aspect_ratio=1., 
            ylabel="Genomic position [Mb]", colorbartitle="Contact frequency",
            yticks=yticks_Mb(N), xticks=xticks_Mb(N))
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
                xticks=xticks_Mb_two_chromosomes(N), yticks=yticks_Mb_two_chromosomes(N),
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
                yticks=yticks_Mb(N), xticks=xticks_Mb(N))
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
            xticks=xticks_Mb_two_chromosomes(N), yticks=yticks_Mb_two_chromosomes(N),
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
            yticks=yticks_Mb(N), xticks=xticks_Mb(N),
            xlabel="Genomic position, Strand 1 [Mb]", ylabel="Genomic position, Strand 2 [Mb]", colorbartitle="Average distance [μm]")

    vline!([N/2+fork[1], fork[2]-N/2], color=:red, linestyle=:dash)
    hline!([N/2+fork[1], fork[2]-N/2], color=:red, linestyle=:dash)
end


"""
Given average long axis positions from simulations with and without smcs,
plot them on the same figure, with an indication of the cell size and fork positions
"""
function plot_compare_av_z(av1, av2, fork, N=404; lbl1="No loop-extruders", lbl2="Loop-extruders", line1=:dash, line2=:solid, plot_forks=true)
        ter=ceil(Int,N/2)

        unrepl1=fork[1]:ter-1
        unrepl2=ter:fork[2]

        repl=vcat(fork[2]:N, 1:fork[1])
        repl_inds= repl.-(ter)
        repl_inds[repl_inds.<=0].+=N+1

        L=fork_to_height(fork)

        old=av2[1:N]
        new=av2[N+1:end]

        #plot the cell height as an outline
        hline([L/2], ribbon=(L, 0), color=:turquoise, fillalpha=0.2)
        hline!([-L/2], color=:turquoise, xticks=xticks_Mb(N))
        hline!([0], color=:grey, linealpha=0.3)

        #plot means and errors for the first set of simulations
        plot!(repl_inds,old[repl], ylabel="Mean long axis position [μm]",xlabel="Genomic position [Mb]", color=1, linestyle=line2)
        plot!(unrepl2.-(ter),old[unrepl2], color=:black,label=lbl2, linestyle=line2)
        plot!(unrepl1.-(ter-N+1),old[unrepl1], color=:black, linestyle=line2)
        plot!(repl_inds,new[repl], color=2, size=[700,400], linestyle=line2)

        old2=av1[1:N]
        new2=av1[N+1:end]

        #plot means and errors for the second set of simulations
        plot!(repl_inds,old2[repl], color=1, linestyle=line1, ylims=(-3.45/2, 3.45/2))
        plot!(unrepl2.-(ter),old2[unrepl2], color=:black,label=lbl1, linestyle=line1)
        plot!(unrepl1.-(ter-N),old2[unrepl1], color=:black,linestyle=line1)
        plot!(repl_inds,new2[repl], color=2, linestyle=line1)
        if plot_forks
                xs=[N/2+fork[1], fork[2]-N/2, N/2+fork[1], fork[2]-N/2]
                ys=[av1[fork[1]], av1[fork[2]], av2[fork[1]+N], av2[fork[2]]]

                scatter!(xs,ys, color=:red, markersize=5, label="Replication fork", markerstrokecolor=:red)
        end
end

"""
Given average long axis positions from simulations with and without smcs,
plot them on the same figure, with an indication of the cell size and fork positions
"""
function plot_compare_av_z_w_error(av1, av2, er1, er2, fork, N=404; lbl1="No loop-extruders", lbl2="Loop-extruders", line1=:dash, line2=:solid, plot_forks=true)
        ter=ceil(Int,N/2)

        unrepl1=fork[1]:ter-1
        unrepl2=ter:fork[2]

        repl=vcat(fork[2]:N, 1:fork[1])
        repl_inds= repl.-(ter)
        repl_inds[repl_inds.<=0].+=N+1

        L=fork_to_height(fork)

        old=av2[1:N]
        new=av2[N+1:end]
        er_old=er2[1:N]
        er_new=er2[N+1:end]

        #plot the cell height as an outline
        hline([L/2], ribbon=(L, 0), color=:turquoise, fillalpha=0.2)
        hline!([-L/2], color=:turquoise, xticks=xticks_Mb(N))
        hline!([0], color=:grey, linealpha=0.3)

        #plot means and errors for the first set of simulations
        plot!(repl_inds,old[repl], ribbon=er_old[repl], ylabel="Mean long axis position [μm]",xlabel="Genomic position [Mb]", color=1, linestyle=line2)
        plot!(unrepl2.-(ter),old[unrepl2], ribbon=er_old[unrepl2],color=:black,label=lbl2, linestyle=line2)
        plot!(unrepl1.-(ter-N+1),old[unrepl1], ribbon=er_old[unrepl1], color=:black, linestyle=line2)
        plot!(repl_inds,new[repl], ribbon=er_new[repl], color=2, size=[700,400], linestyle=line2)

        old2=av1[1:N]
        new2=av1[N+1:end]
        er_old=er1[1:N]
        er_new=er1[N+1:end]

        #plot means and errors for the second set of simulations
        plot!(repl_inds,old2[repl], ribbon=er_old[repl], color=1, linestyle=line1, ylims=(-3.45/2, 3.45/2))
        plot!(unrepl2.-(ter),old2[unrepl2], ribbon=er_old[unrepl2], color=:black,label=lbl1, linestyle=line1)
        plot!(unrepl1.-(ter-N),old2[unrepl1],ribbon=er_old[unrepl1], color=:black,linestyle=line1)
        plot!(repl_inds,new2[repl], ribbon=er_new[repl], color=2, linestyle=line1)
        if plot_forks
                xs=[N/2+fork[1], fork[2]-N/2, N/2+fork[1], fork[2]-N/2]
                ys=[av1[fork[1]], av1[fork[2]], av2[fork[1]+N], av2[fork[2]]]

                scatter!(xs,ys, color=:red, markersize=5, label="Replication fork", markerstrokecolor=:red)
        end
end

"""
Given average long axis positions from simulations at different time points,
plot them on the same figure
"""
function plot_av_zs_all(avs, forks, N=405)
        ter=ceil(Int,N/2)
        nt=size(avs)[2]
        if N>1000
                mon_kb=1
        else
                mon_kb=10
        end
        hline([0], color=:grey, linealpha=0.3, ylims=(-1.5, 1.5), xticks=xticks_Mb(N))

        for i in 1:nt
                α=1-(i-1)/(nt)

                unrepl1=forks[1,i]:ter-1
                unrepl2=ter:min(N,forks[2,i])

                repl=vcat(min(N,forks[2,i]-1):N, 1:forks[1,i]+1)
                repl_inds= repl.-(ter)
                repl_inds[repl_inds.<=0].+=N+1

                old=avs[1:N,i]
                new=avs[N+1:end,i]
                t=Int(round(fork_to_time(forks[:,i],N=N, mon_kb=mon_kb), digits=-1))
                if i==1
                        unrepl2=vcat(ter:N,1:ter)
                        plot!(old[unrepl2], color=:black, label="0 min")
                else 
                        plot!(repl_inds,old[repl], ylabel="Mean long axis position [μm]",xlabel="Genomic position [Mb]", color=1, alpha=α)
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
        plot!(zs[1,:], ribbon=stds[1,:], label="Ori 1", color=1, ylabel="Mean long axis position [μm]")
        plot!(zs[N+1,:], ribbon=stds[N+1,:], label="Ori 2", color=2, xlabel="Replicated length [Mb]")
        plot!(zs[ter,:], ribbon=stds[ter,:], label="Ter", color=:black, alpha=0.8, xticks=xticks_Mb(N, false))
end

"""
Plot a kymograph of the long axis positions of the centers of mass of the chromosomal strands as a function of time
Include an indicator for the cell size as a function of time
"""
function com_kymograph(zs, stds)
        N=size(zs)[2]-1
        hs=map(i->R_to_height(i), 0:N)

        plot(hs./2, ribbon=(hs, zeros(N+1)), color=:turquoise, fillalpha=0.2)
        plot!(zs[1,:], ribbon=stds[1,:], label="Strand 1", color=1, ylabel="Mean center of mass\nlong axis position [μm]")
        plot!(zs[2,:], ribbon=stds[2,:], label="Strand 2", color=2, xlabel="Replicated length [Mb]")
        plot!(zs[3,:], ribbon=stds[3,:], label="Ter strand", color=:black, alpha=0.8, xticks=xticks_Mb(N, false))
        plot!(-hs./2, color=:turquoise)
end

"""
Given a set of replicated lengths and segregated fractions, plot a scatter of the fractions
with error bars

segregated_fractions should be an array of size (2, n) where n is the number of replicated lengths
The first row should be the segregated fractions without loop-extruders, the second row with loop-extruders
"""
function plot_compare_segregated_fractions(Rs, segregated_fractions, segregated_fractions_std; lbl2="Loop-extruders", lbl1="No loop-extruders", lbl3="Ideal chain")
        Rs_prepend = vcat(0, Rs)
        
        segregated_fractions_prepend = vcat(0, segregated_fractions)
    
        # First, draw the lines without any markers
        p = plot(Rs_prepend, segregated_fractions_prepend, xticks=xticks_Mb(404, false),
                 ylabel="Mean segregated fraction", xlabel="Replicated length [Mb]",
                 size=(800,400), ylims=(0,1), xlims=(0,415), color=[:teal :black :magenta], 
                 fillalpha=[0.2 0.4 0.4])
    
        # Now, add the markers on top but skip the first one
        scatter!(Rs, segregated_fractions, yerror=segregated_fractions_std,
                 marker=[:circle :square :diamond], markersize=10, color=[:teal :black :magenta],
                 label=[lbl1 lbl2 lbl3], msc=color=[:teal :black :magenta])    
        return p
end
    


"""
Given long vectors of segregated fractions and their standard deviations as a function of R from dynamic simulations
Plot the segregated fractions as a function of R

segregated_fractions_1 should correspond to data without loop-extruders, segregated_fractions_2 with loop-extruders
"""
function plot_compare_segregated_fractions(segregated_fractions, segregated_fractions_std; N=404, std_ideal=false, lbl1="No loop-extruders", lbl2="Loop-extruders", lbl3="Ideal chain")
        ribbons=copy(segregated_fractions_std)
        if !std_ideal
                ribbons[:,3].=NaN
        end
        plot(1:N, segregated_fractions, ribbon=ribbons, ylims=(0,1), 
                xlabel="Replicated length [Mb]", ylabel="Mean segregated fraction",
                xticks=xticks_Mb(N,false), label=[lbl1 lbl2 lbl3],  
                color=[:teal :black :magenta], size=(800,400), fillalpha=[0.4 0.2 0.4])
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
        plot!(long_axis_positions[1,:], xlabel="Simulation time", ylabel="Mean long axis position [μm]", label="Ori 1", color=1)
        plot!(long_axis_positions[N+1,:], label="Ori 2", ylims=(-L_max/2,L_max/2), legend=:outerright, color=2)
        plot!(long_axis_positions[Int(N/2),:], label="Ter", xticks=:auto, color=:black)
end

"""
Plot that helps check where SMCs are.
"""
function smc_distr(smcs, R, N)
        scatter(smcs[1,:], smcs[2,:], xlims=(0,2*N), ylims=(0,2*N))
        vline!([N,2*N], color=:black)
        hline!([N,2*N], color=:black)
        vline!([R/2, N+R/2], color=:red, linestyle=:dash)
        hline!([R/2, N+R/2], color=:red, linestyle=:dash)
end