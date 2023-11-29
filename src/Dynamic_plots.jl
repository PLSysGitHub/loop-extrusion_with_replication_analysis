include("Plot_functions.jl")

################################################################################
# Pipeline for loading data and making plots.
# To adjust the actual plots, edit Plot_functions.jl
################################################################################

"""
Create plots of Hi-C maps and distance maps for given time points in
dynamic simulations
"""
function make_time_point_plots(times, data_dir,plot_dir, monomer_size; N=404, filetype="png")
        for subf in ["Hi_C/all", "Hi_C/strand_1", "Hi_C/strand_2", "Hi_C/trans", "Chipseq", "Distance_matrix", "Distance_trans", "Mean_z"]
                if !isdir(plot_dir*subf)
                        mkpath(plot_dir*subf)
                end
        end

        if filetype[1]!='.'
                filetype="."*filetype
        end

        subdirs=readdir(data_dir, join=true)
        subdirs=subdirs[isdir.(subdirs)]

        for subdir in subdirs
                files=readdir(subdir, join=true)
                deleteat!(files, contains.(files, "segregation.h5"))
                for filename in files
                        #each figure saved in two places; folder for single simulation
                        out_dir=replace(filename, data_dir=>plot_dir)
                        out_dir=replace(out_dir, ".h5"=>"/")

                        println("Making plots for $filename...")

                        #other folders by plot type
                        out_ext=replace(filename, data_dir=>"")
                        out_ext=replace(out_ext, "/"=>"_")
                        out_ext=replace(out_ext, ".h5"=>"")

                        h5open(filename, "r") do file
                                forks=floor.(Int,read(file, "mean_forks"))

                                #individual hic and chipseq maps for each time point
                                if !contains(subdir, "No_smcs") && !contains(subdir, "nosmcs") && !contains(subdir, "Ideal")
                                        chipseq=read(file, "mean_chipseq")
                                end
                                contacts=read(file, "mean_hic")
                                distance_maps=read(file, "av_dist")*monomer_size
                                mean_zs=read(file, "mean_zs")*monomer_size
                                plot_av_zs_all(mean_zs, forks, N)
                                savefig(plot_dir*"Mean_z/"*out_ext*filetype)
                                for (index,t) in enumerate(times)
                                        file_prefix_shared_dir=out_ext*"_t_$(t)_"
                                        R=forks[1,index]+N-forks[2,index]
                                        #Hi-C maps

                                        old=contacts[:,:,2, index]
                                        new=contacts[:,:,3, index]
                                        inter=contacts[:,:,4, index]

                                        old.+=transpose(old)
                                        new.+=transpose(new)

                                        plot_hic(old, new, inter, R=R)
                                        savefig(plot_dir*"Hi_C/all/"*file_prefix_shared_dir*filetype)

                                        plot_hic(old, R=R)
                                        savefig(plot_dir*"Hi_C/strand_2/"*file_prefix_shared_dir*filetype)

                                        plot_hic(new, R=R)
                                        savefig(plot_dir*"Hi_C/strand_1/"*file_prefix_shared_dir*filetype)

                                        plot_hic(inter, R=R)
                                        savefig(plot_dir*"Hi_C/trans/"*file_prefix_shared_dir*filetype)

                                        #Distance maps
                                        plot_distance_matrix(distance_maps[:,:,index], forks[:,index])
                                        savefig(plot_dir*"Distance_matrix/"*file_prefix_shared_dir*filetype)

                                        plot_distance_matrix_oldnew(distance_maps[:,:,index], forks[:,index])
                                        savefig(plot_dir*"Distance_trans/"*file_prefix_shared_dir*filetype)

                                        #Chip-seq
                                        if !contains(subdir, "No_smcs") && !contains(subdir, "nosmcs") && !contains(subdir, "Ideal")
                                                plot(chipseq[:,index], xlabel="Genomic position [Mb]", ylabel="Mean number loop-extruders", 
                                                        color=:black, ylims=(0,1.1), size=[500,350], xticks=(0:200:400,0:2:4))
                                                vline!(forks[:,index], color=:red, linestyle=:dash)
                                                savefig(plot_dir*"Chipseq/"*file_prefix_shared_dir*filetype)
                                        end
                                end

                        end
                end
        end
end

"""
Create plots of segregated fraction over time as well as ori kymographs for dynamic simulations.
"""
function make_time_course_plots(data_dir,plot_dir, monomer_size, N; filetype="png")
        for subf in ["Segregated_fractions", "Ori_kymograph"]
                if !isdir(plot_dir*subf)
                        mkpath(plot_dir*subf)
                end
        end

        if filetype[1]!='.'
                filetype="."*filetype
        end

        subdirs=readdir(data_dir, join=true)
        subdirs=subdirs[isdir.(subdirs)]

        for subdir in subdirs
                files=readdir(subdir, join=true)
                files=files[contains.(files,"segregation")]
                for filename in files

                        println("Making plots for $filename...")

                        #other folders by plot type
                        out_ext=replace(filename, data_dir=>"")
                        out_ext=replace(out_ext, "/"=>"_")
                        out_ext=replace(out_ext, ".h5"=>"")

                        h5open(filename, "r") do file
                                zs=read(file,"zs")*monomer_size
                                zs_std=read(file, "zs_std")*monomer_size
                                ori_ter_kymograph(zs,zs_std)
                                savefig(plot_dir*"Ori_kymograph/"*out_ext*filetype)

                                seg_inds=read(file, "seg_index")[10:end]
                                seg_std=read(file, "seg_index_std")[10:end]
                                
                                plot(seg_inds, ribbon=seg_std, ylims=(0,1), xlabel="Replicated length [Mb]", color=:black, fillalpha=0.2,
                                        ylabel="Segregated fraction", xticks=(1:100:404, string.([0,1, 2, 3, 4])))
                                savefig(plot_dir*"Segregated_fractions/"*out_ext*filetype)
                        end
                end
        end
end

"""
Compare segregated fraction and radial distance over time for simulations with and without smcs.
"""
function plot_compare_segregation_data(file_no_smc, file_smc, file_ideal, plot_dir, monomer_size, N; filetype="png", lbl1="No loop-extruders", lbl2="Loop-extruders", lbl3="Ideal chain")
        if !isdir(plot_dir)
                mkpath(plot_dir)
        end

        if filetype[1]!='.'
                filetype="."*filetype
        end
        
        seg_inds=zeros(N,3)
        seg_std=zeros(N,3)

        for (index, file) in enumerate([file_no_smc, file_smc, file_ideal])
                if isfile(file)
                        h5open(file) do f
                                seg_inds[:,index].=read(f, "seg_index")[1:N]
                                seg_std[:,index].=read(f, "seg_index_std")[1:N]
                        end
                else
                        seg_inds[:,index].=NaN
                        seg_std[:,index].=NaN
                end
        end

        plot_compare_segregated_fractions(seg_inds, seg_std, lbl1=lbl1, lbl2=lbl2, lbl3=lbl3)
        savefig(plot_dir*"compare_segregated_fractions"*filetype)

        plot_compare_segregated_fractions(seg_inds, seg_std, std_ideal=true, lbl1=lbl1, lbl2=lbl2, lbl3=lbl3)
        savefig(plot_dir*"compare_segregated_fractions_ribbon_ideal"*filetype)
end

"""
Compare mean long axis position of monomers at given time points for simulations with and without smcs.
"""
function plot_compare_time_point_data(times, file_no_smc, file_smc, plot_dir, monomer_size, N;skip_done=true, filetype="png", lbl1="No loop-extruders", lbl2="Loop-extruders", lbl3="Ideal chain", line1=:dash, line2=:solid)
        if !isdir(plot_dir)
                mkpath(plot_dir)
        end
        if filetype[1]!='.'
                filetype="."*filetype
        end

        h5open(file_no_smc) do f1
                forks=floor.(Int,read(f1, "mean_forks"))
                zs_no_smc=read(f1, "mean_zs")*monomer_size

                h5open(file_smc) do f2
                        zs_smc=read(f2, "mean_zs")*monomer_size

                        for (index,t) in enumerate(times)
                                if t>0
                                        plot_compare_av_z(zs_no_smc[:,index], zs_smc[:,index], forks[:,index], lbl1=lbl1, lbl2=lbl2,line1=line1, line2=line2)
                                        savefig(plot_dir*"compare_average_long_axis_pos_t_$t"*filetype)
                                end
                        end
                end
        end

end