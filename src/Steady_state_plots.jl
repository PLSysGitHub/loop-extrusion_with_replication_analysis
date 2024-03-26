include("Plot_functions.jl")

"""
Plot Hi-C maps and distance maps for steady state simulations
"""
function make_slow_eq_plots(data_dir,plot_dir, monomer_size; N=404, filetype="png", skip_done=true)
    for subf in ["Hi_C/all", "Hi_C/strand_1", "Hi_C/strand_2", "Hi_C/trans", "Distance_matrix", "Distance_trans"]
        if !isdir(plot_dir*subf)
            mkpath(plot_dir*subf)
        end
    end

    if filetype[1]!='.'
        filetype="."*filetype
    end

    for subdir in subdirs(data_dir)
        for subsubdir in subdirs(subdir)
            filename=subsubdir*"/slow_steady_state_stats.h5"
                        
            println("Making plots for $filename...")

            #folders by plot type; filename for single simulation
            out_ext=replace(filename, data_dir=>"")
            out_ext=replace(out_ext, "/"=>"_")
            out_ext=replace(out_ext, ".h5"=>"")

            if isfile(plot_dir*"Distance_trans/"*out_ext*filetype) && skip_done
                println("Found plot $filename, skipping...")
                continue
            end

            h5open(filename, "r") do file
                R=parse_after(filename)
                fork=floor.(Int,read(file, "fork"))

                contacts=read(file, "mean_hic")
                distance_maps=read(file, "av_dist")*monomer_size

                #Hi-C maps
                old=contacts[:,:,2]
                new=contacts[:,:,3]
                inter=contacts[:,:,4]

                old.+=transpose(old)
                new.+=transpose(new)

                pl=plot_hic(old, new, inter, R=R)
                savefig(pl, plot_dir*"Hi_C/all/"*out_ext*filetype)

                pl=plot_hic(old, R=R)
                savefig(pl, plot_dir*"Hi_C/strand_1/"*out_ext*filetype)

                pl=plot_hic(new, R=R)
                savefig(pl, plot_dir*"Hi_C/strand_2/"*out_ext*filetype)

                pl=plot_hic(inter, R=R)
                savefig(pl, plot_dir*"Hi_C/trans/"*out_ext*filetype)

                #Distance maps
                pl=plot_distance_matrix(distance_maps, fork)
                savefig(pl, plot_dir*"Distance_matrix/"*out_ext*filetype)

                pl=plot_distance_matrix_oldnew(distance_maps, fork)
                savefig(pl, plot_dir*"Distance_trans/"*out_ext*filetype)
            end
        end
    end
end

"""
Plot segregated fraction and long axis position of oris over time for steady state simulations.

Plots can be used to assess whether simulations have converged.
"""
function make_convergence_plots(data_dir, plot_dir, monomer_size; N=404, spacer=1, filetype="png")
    if !isdir(plot_dir)
        mkpath(plot_dir)
    end

    if filetype[1]!='.'
        filetype="."*filetype
    end

    for dir in subdirs(data_dir)
        files=readdir(dir, join=true)
        files=files[contains.(files, "z_over_sims.txt")]
        out_plot=replace(dir, data_dir=>plot_dir)
        if !ispath(out_plot)
            mkpath(out_plot)
        end
        for f in files
            segr_file=replace(f, "z_over_sims"=>"seg_ind_over_sims")
            R=parse_after(f)

            if R>10
                zs=readdlm(f)*monomer_size
                seg_inds=vec(readdlm(segr_file))

                last_index= any(isnan.(seg_inds)) ? findfirst(isnan,seg_inds)-1 : length(seg_inds)

                plot_convergence_oris(zs[:,1:last_index],R)
                savefig(replace(replace(f, dir=>out_plot), ".txt"=>"_ori_positions")*filetype)

                plot(seg_inds[1:last_index], xlabel="Simulation time", ylabel="Mean segregated fraction", ylims=(0,1), xticks=:auto, color=:black)
                savefig(replace(replace(f, dir=>out_plot), ".txt"=>"_segregated_fractions")*filetype)
            end
        end
    end
end

"""
Compare steady state data for simulations with and without SMCs, at given values of the replicated length R.
"""
function compare_point_data(Rs::Vector{Int}, dir_no_smc::String, dir_smc::String, dir_ideal::String, plot_dir::String, 
        monomer_size::Number, N::Int; filetype="png", lbl1="No loop-extruders", lbl2="Loop-extruders", lbl3="Ideal chain", 
        line1=:dash, line2=:solid, parSval=4040.0, P_U=0.0)
    if !isdir(plot_dir)
            mkpath(plot_dir)
    end

    if filetype[1]!='.'
        filetype="."*filetype
    end

    files_no_smc=readdir(dir_no_smc, join=true)
    files_no_smc=files_no_smc[isdir.(files_no_smc)].*"/steady_state_stats.h5"

    files_smc=readdir(dir_smc, join=true)
    files_smc=files_smc[isdir.(files_smc)].*"/steady_state_stats.h5"

    #Added new data with different loading rates and P_U; make sure loading the right ones
    if contains(files_smc[1], "parSstrength_")
        files_smc=filter(x->contains(x, "parSstrength_$(parSval)_"), files_smc)
    end
    if contains(files_smc[1], "stallFork")
        files_smc=filter(x->contains(x, "stallFork_$(P_U)_"), files_smc)
    end

    files_ideal=readdir(dir_ideal, join=true)
    files_ideal=files_ideal[isdir.(files_ideal)].*"/steady_state_stats.h5"
    
    segregated_fractions=zeros(length(Rs),3)
    segregated_fractions_std=zeros(length(Rs),3)

    # Loop over R and read segregated fractions. Plot comparison of long axis positions.
    for (index,R) in enumerate(Rs)
        file_no_smc=files_no_smc[findfirst(contains("R_$(R)_"), files_no_smc)]
        file_smc=files_smc[findfirst(contains("R_$(R)_"), files_smc)]
        file_ideal=""
        try
            file_ideal=files_ideal[findfirst(contains("R_$(R)_"), files_ideal)]
        catch
            println("No ideal polymer file for R=$R")
        end
        fork=[Int(R/2), Int(N-R/2)]

        h5open(file_no_smc) do f1
                zs_no_smc=read(f1, "mean_zs")*monomer_size
                num_samp_no_smc=read(f1, "num_samp")
                zs_no_smc_error=read(f1, "zs_std")*monomer_size/sqrt(num_samp_no_smc)

                segregated_fractions[index,1]=read(f1, "seg_index")
                segregated_fractions_std[index,1]=read(f1, "seg_index_std")

                h5open(file_smc) do f2
                    zs_smc=read(f2, "mean_zs")*monomer_size
                    num_samp_smc=read(f2, "num_samp")
                    zs_smc_error=read(f2, "zs_std")*monomer_size/sqrt(num_samp_smc)

                    segregated_fractions[index,2]=read(f2, "seg_index")
                    segregated_fractions_std[index,2]=read(f2, "seg_index_std")

                    plot_compare_av_z(zs_no_smc, zs_smc, fork, lbl1=lbl1, lbl2=lbl2, line1=line1, line2=line2)
                    savefig(plot_dir*"compare_average_long_axis_positions_R_$R"*filetype)

                    plot_compare_av_z_w_error(zs_no_smc, zs_smc,zs_no_smc_error, zs_smc_error, fork, lbl1=lbl1, lbl2=lbl2, line1=line1, line2=line2, plot_forks=false)
                    savefig(plot_dir*"with_error_compare_average_long_axis_positions_R_$R"*filetype)
                end
        end

        if file_ideal!=""
            h5open(file_ideal) do f3
                segregated_fractions[index,3]=read(f3, "seg_index")
                segregated_fractions_std[index,3]=read(f3, "seg_index_std")
            end
        else
            segregated_fractions[index,3]=NaN
            segregated_fractions_std[index,3]=NaN
        end
    end

    plot_compare_segregated_fractions(Rs, segregated_fractions, segregated_fractions_std, lbl1=lbl1, lbl2=lbl2, lbl3=lbl3)
    savefig(plot_dir*"compare_segregated_fractions"*filetype)
end

"""
Read the chipseq profiles and make plots
"""
function make_chipseq_plots(data_dir, plot_dir;N=404, filetype="png")
    if filetype[1]!='.'
        filetype="."*filetype
    end
    for dir in subdirs(data_dir)
        if contains(dir, "No_smcs") || contains(dir, "Ideal")
            continue
        end
        if !isdir(replace(dir, data_dir=>plot_dir))
            mkpath(replace(dir, data_dir=>plot_dir))
        end
        for subdir in subdirs(dir)
            filename=subdir*"/steady_state_chipseq.h5"
            out_plot=replace(subdir, data_dir=>plot_dir)
            h5open(filename, "r") do file
                R=parse_after(filename)
                chipseq=read(file, "mean_chipseq")
                forks=read(file, "fork")

                plot(chipseq, xlabel="Genomic position [Mb]", ylabel="Mean number loop-extruders", 
                    color=:black, ylims=(0,1.1), size=[500,350], xticks=(0:200:400,0:2:4))
                vline!(forks, color=:red, linestyle=:dash)

                savefig(out_plot*"chipseq_profile"*filetype)

                smc_positions=read(file, "mean_smc_positions")
                heatmap(smc_positions, xlabel="Genomic position [Mb]", ylabel="Genomic position [Mb]", 
                    color=:viridis, size=[500,350], xticks=(0:200:800), yticks=(0:200:800), clims=(0,0.002))
                vline!(forks, color=:red, linestyle=:dash)
                hline!(forks, color=:red, linestyle=:dash)
                savefig(out_plot*"smc_positions_heatmap"*filetype)
            end
        end
    end
end