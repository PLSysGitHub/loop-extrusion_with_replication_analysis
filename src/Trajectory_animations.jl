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

function animate_file(in_file_name::String, times, out_file_name::String; N=404, spacer=1, monomer_size=0.129)
    nt=length(times)
    h5open(in_file_name) do file
        list_inds=get_indices_at_times(file,times)
        anim = @animate for time_ind in 1:nt
            pos=fetch_raw_pos_at_ind(file, list_inds[time_ind], spacer=spacer, N=N)*monomer_size
            forks=fetch_fork_at_ind(file, list_inds[time_ind],spacer=spacer)

            plot_strands(pos,forks,times[time_ind])
        end
        gif(anim, out_file_name, fps=10)
    end
end

function snapshots_file(in_file_name::String, times, out_file_name::String; N=404, spacer=1, monomer_size=0.129, file_type=".pdf")
    nt=length(times)
    if file_type[1]!='.'
        file_type="."*file_type
    end
    h5open(in_file_name) do file
        list_inds=get_indices_at_times(file,times)
        for time_ind in 1:nt
            pos=fetch_raw_pos_at_ind(file, list_inds[time_ind], spacer=spacer, N=N)*monomer_size
            forks=fetch_fork_at_ind(file, list_inds[time_ind],spacer=spacer)

            plot_strands(pos,forks,times[time_ind])
            savefig(out_file_name*"_$time_ind"*file_type)
        end
    end
end

function one_animation_per_sim(data_parent_dir="Simulation_data/Dynamic/", out_dir="Animations/"; times=0:0.5:80, which=80, skip_done=true, N=N, spacer=1, monomer_size=0.129)
    if !isdir(out_dir)
        mkpath(out_dir)
    end
    subdirectories=subdirs(data_parent_dir)
    for subdir in subdirectories #No_smcs, Topological etc
        dirs=subdirs(subdir)
        for dir in dirs
            out_file_name=replace(dir, data_parent_dir=>"")
            out_file_name=out_dir*replace(out_file_name, "/"=>"_")*"_$which.gif"
            if skip_done && isfile(out_file_name)
                println("Found file $out_file_name. Skipping animation.")
            else
                if length(readdir(dir)) >= which
                    filename=readdir(dir, join=true)[which]
                    animate_file(filename, times, out_file_name, N=N, spacer=spacer, monomer_size=monomer_size)
                else
                    println("$dir has only $(length(readdir(dir))), no animation for $which...")
                end
            end
        end
    end
end

function snapshots_per_sim(data_parent_dir="Simulation_data/Dynamic/", out_parent_dir="Snapshots/"; skip_done=true, times=0:5:80, which=80, N=404, file_type=".pdf", spacer=1, monomer_size=0.129)

    subds=subdirs(data_parent_dir) #All dynamic simulations
    for subdir in subds #No_smcs, Topological etc
        dirs=subdirs(subdir)
        for dir in dirs #Different variables
            out_dir=out_parent_dir*replace(dir, data_parent_dir=>"")
            if !isdir(out_dir)
                mkpath(out_dir)
            end
            out_file_name=out_dir*"/file_$(which)_snapshot"

            if isfile(out_file_name) && skip_done #already done
                continue
            else
                if length(readdir(dir)) >= which #there are enough files to make snapshots for requested run
                    filename=readdir(dir, join=true)[which]
                    snapshots_file(filename, times, out_file_name, N=N, spacer=spacer, monomer_size=monomer_size, file_type=file_type)
                else
                    println("$dir has only $(length(readdir(dir))), no snapshots for $which...")
                end
            end
        end
    end
end