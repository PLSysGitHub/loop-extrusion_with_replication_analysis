include("Plot_functions.jl")

"""
Plot Hi-C maps and distance maps for steady state simulations
"""
function make_slow_eq_plots(data_dir,plot_dir, monomer_size; N=404)
    for subf in ["Hi_C/all", "Hi_C/strand_1", "Hi_C/strand_2", "Hi_C/trans", "Distance_matrix", "Distance_trans"]
        if !isdir(plot_dir*subf)
            mkpath(plot_dir*subf)
        end
    end

    for subdir in subdirs(data_dir)
        for subsubdir in subdirs(subdir)
            filename=subsubdir*"/slow_steady_state_stats.h5"
                        
            println("Making plots for $filename...")

            #folders by plot type; filename for single simulation
            out_ext=replace(filename, data_dir=>"")
            out_ext=replace(out_ext, "/"=>"_")
            out_ext=replace(out_ext, ".h5"=>"")

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
                png(pl, plot_dir*"Hi_C/all/"*out_ext)

                pl=plot_hic(old, R=R)
                png(pl, plot_dir*"Hi_C/strand_1/"*out_ext)

                pl=plot_hic(new, R=R)
                png(pl, plot_dir*"Hi_C/strand_2/"*out_ext)

                pl=plot_hic(inter, R=R)
                png(pl, plot_dir*"Hi_C/trans/"*out_ext)

                #Distance maps
                pl=plot_distance_matrix(distance_maps, fork)
                png(pl, plot_dir*"Distance_matrix/"*out_ext)

                pl=plot_distance_matrix_oldnew(distance_maps, fork)
                png(pl, plot_dir*"Distance_trans/"*out_ext)
            end
        end
    end
end

"""
Plot segregated fraction and long axis position of oris over time for steady state simulations.

Plots can be used to assess whether simulations have converged.
"""
function make_convergence_plots(data_dir, plot_dir, monomer_size; N=404, spacer=1)
    if !isdir(plot_dir)
        mkpath(plot_dir)
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
                png(replace(replace(f, dir=>out_plot), ".txt"=>"_ori_positions"))

                plot(seg_inds[1:last_index], xlabel="Simulation time", ylabel="Segregated fraction", ylims=(0,1), xticks=:auto, color=:black)
                png(replace(replace(f, dir=>out_plot), ".txt"=>"_segregated_fractions"))
            end
        end
    end
end

"""
Compare steady state data for simulations with and without SMCs, at given values of the replicated length R.
"""
function compare_point_data(Rs::Vector{Int}, dir_no_smc::String, dir_smc::String, plot_dir::String, monomer_size::Number, N::Int)
    if !isdir(plot_dir)
            mkpath(plot_dir)
    end

    files_no_smc=readdir(dir_no_smc, join=true)
    files_no_smc=files_no_smc[isdir.(files_no_smc)].*"/steady_state_stats.h5"

    files_smc=readdir(dir_smc, join=true)
    files_smc=files_smc[isdir.(files_smc)].*"/steady_state_stats.h5"
    
    segregated_fractions=zeros(2,length(Rs))
    segregated_fractions_std=zeros(2,length(Rs))

    # Loop over R and read segregated fractions. Plot comparison of long axis positions.
    for (index,R) in enumerate(Rs)
        file_no_smc=files_no_smc[findfirst(contains("R_$(R)_"), files_no_smc)]
        file_smc=files_smc[findfirst(contains("R_$(R)_"), files_smc)]

        fork=[Int(R/2), Int(N-R/2)]

        h5open(file_no_smc) do f1
                zs_no_smc=read(f1, "mean_zs")*monomer_size

                segregated_fractions[1,index]=read(f1, "seg_index")
                segregated_fractions_std[1,index]=read(f1, "seg_index_std")

                h5open(file_smc) do f2
                    zs_smc=read(f2, "mean_zs")*monomer_size

                    segregated_fractions[2,index]=read(f2, "seg_index")
                    segregated_fractions_std[2,index]=read(f2, "seg_index_std")

                    plot_compare_av_z(zs_no_smc, zs_smc, fork)
                    png(plot_dir*"compare_average_long_axis_positions_R_$R")
                end
        end
    end

    plot_compare_segregated_fractions(Rs, segregated_fractions, segregated_fractions_std)
    png(plot_dir*"compare_segregated_fractions")
end