using Plots
include("helpers.jl")
using Plots.PlotMeasures
gr(palette=:tab10)

function plot_strands(pos, forks, t)
    h=3.3 #maximum height in μm
    x=pos[1,:]
    y=pos[2,:]
    z=pos[3,:]

    N=Int(length(x)/2)
    old=vcat(round(Int,forks[2]):N, 1:round(Int,forks[1]))
    new=old.+N
    unreplicated=max(1,round(Int,forks[1])):min(round(Int,forks[2]), N)

    plot3d(z[old],x[old],y[old], color=1,label="Old", title="$t min",
        xlims=(-h/2,h/2), ylims=(-h/2,h/2),zlims=(-h/2,h/2), xlabel="", 
        framestyle=:none, ticks=false, grid=:hide, showaxis=:hide, 
        legend=false, linewidth=1.6, size=(600,600), margin=0mm)
    plot3d!(z[new],x[new],y[new], color=2, label="New")
    plot3d!(z[unreplicated],x[unreplicated],y[unreplicated], color=:black, label="Unreplicated")

    return current()
end

function animate_file(in_file_name::String, times, out_file_name::String, N=404)
    nt=length(times)
    h5open(in_file_name) do file
        if N==404
            spacer=1
            monomer_size=0.129
        else
            spacer=Int(N/405)
            N=405
            monomer_size=0.088
        end
        list_inds=get_indices_at_times(file,times)
        anim = @animate for time_ind in 1:nt
            pos=fetch_raw_pos_at_ind(file, list_inds[time_ind], spacer=spacer, N=N)*monomer_size
            forks=fetch_fork_at_ind(file, list_inds[time_ind],spacer=spacer)

            # if subdir!=data_parent_dir*"No_smcs"
            #     SMCs=fetch_smcs_at_inds(file, ind)
            # end

            plot_strands(pos,forks,times[time_ind])
        end
        gif(anim, out_file_name, fps=10)
    end
end

function snapshots_file(in_file_name::String, times, out_file_name::String, N=404)
    nt=length(times)
    h5open(in_file_name) do file
        if N==404
            spacer=1
            monomer_size=0.129
        else
            spacer=Int(N/405)
            N=405
            monomer_size=0.088
        end
        list_inds=get_indices_at_times(file,times)
        for time_ind in 1:nt
            pos=fetch_raw_pos_at_ind(file, list_inds[time_ind], spacer=spacer, N=N)*monomer_size
            forks=fetch_fork_at_ind(file, list_inds[time_ind],spacer=spacer)

            # if subdir!=data_parent_dir*"No_smcs"
            #     SMCs=fetch_smcs_at_inds(file, ind)
            # end

            plot_strands(pos,forks,times[time_ind])
            png(out_file_name*"_$time_ind")
        end
    end
end

function one_animation_per_sim(data_parent_dir="Simulation_data/Dynamic/", out_dir="Animations/"; times=0:0.5:80, which=80, skip_done=true)
    if !isdir(out_dir)
        mkpath(out_dir)
    end
    subdirs=readdir(data_parent_dir, join=true)
    subdirs=subdirs[isdir.(subdirs)]
    for subdir in subdirs #No_smcs, Topological etc
        dirs=readdir(subdir, join=true)
        dirs=dirs[isdir.(dirs)]
        for dir in dirs
            out_file_name=replace(dir, data_parent_dir=>"")
            out_file_name=out_dir*replace(out_file_name, "/"=>"_")*"_$which.gif"
            if skip_done && isfile(out_file_name)
                println("Found file $out_file_name. Skipping animation.")
            else
                if length(readdir(dir)) >= which
                    filename=readdir(dir, join=true)[which]
                    animate_file(filename, times, out_file_name)
                else
                    println("$dir has only $(length(readdir(dir))), no animation for $which...")
                end
            end
        end
    end
end

function snapshots_per_sim(data_parent_dir="Simulation_data/Dynamic/", out_parent_dir="Snapshots/"; times=0:5:80, which=80)

    subdirs=readdir(data_parent_dir, join=true)
    subdirs=subdirs[isdir.(subdirs)]
    for subdir in subdirs #No_smcs, Topological etc
        dirs=readdir(subdir, join=true)
        dirs=dirs[isdir.(dirs)]
        for dir in dirs
            out_dir=out_parent_dir*replace(dir, data_parent_dir=>"")
            if !isdir(out_dir)
                mkpath(out_dir)
            end
            out_file_name=out_dir*"/file_$(which)_snapshot"

            if length(readdir(dir)) >= which
                filename=readdir(dir, join=true)[which]
                snapshots_file(filename, times, out_file_name)
            else
                println("$dir has only $(length(readdir(dir))), no snapshots for $which...")
            end
        end
    end
end