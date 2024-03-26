"""
This file runs a pipeline that calculates statistics given a folder of loop-extrusion simulation results

Takes a path to input data as an argument.

Calculates:

- Mean z-positions
- Mean fork positions
- Mean local extensions
- Hi-C maps
- Mean chipseq [if there's smcs in the sims]

averaged over all samples, at roughly the times requested in 'times'

The statistics are saved in .h5 files in the directory Stats
"""

include("Statistics_functions.jl")

"""
Calculate stats corresponding to a given replication stages

"""
function calc_time_point_stats(times, N=404,data_parent_dir="Simulation_data/Dynamic/",out_parent_dir="Stats/Dynamic/"; skip_done=true, spacer=1)

    subdirs=readdir(data_parent_dir, join=true)
    subdirs=subdirs[isdir.(subdirs)]

    for subdir in subdirs
        out_dir=replace(subdir, data_parent_dir=>out_parent_dir)
        if !isdir(out_dir)
            mkpath(out_dir)
        end

        simdirs=readdir(subdir, join=true)
        simdirs=simdirs[isdir.(simdirs)]
        for dir in simdirs
            out_file_name=replace(dir, data_parent_dir=>out_parent_dir)*".h5"
            if skip_done && isfile(out_file_name)
                println("Found previous results for $dir, skipping...")
            else
                println("Calculating statistics for $dir")
                num_samp=0
                ter_distinct=true #terminus distinct during dynamic sims
                no_orient=!(contains(dir, "no_tether")) #don't relabel oris if using a tether
                #Initialize output arrays
                mean_zs=zeros(2*N, length(times))
                mean_forks=zeros(2, length(times))
                mean_hic=zeros(N, N, 4, length(times)) #the 4 is the dimension for all/old/new/inter
                av_dist=zeros(2*N,2*N, length(times))
                if !contains(subdir, "No_smcs") && !contains(subdir, "nosmcs")
                    mean_chipseq=zeros(N, length(times))
                end
                #Loop over files, each containing a single time trajectory
                for filename in readdir(dir,join=true)
                        h5open(filename, "r") do file
                            list_inds=get_indices_at_times(file,times)
                            positions, forks=fetch_pos_fork_at_inds(file, list_inds, no_turn=no_orient, N=N, spacer=spacer, ter_distinct=ter_distinct)
                            if !contains(subdir, "No_smcs") && !contains(subdir, "nosmcs") && !contains(subdir, "Ideal")
                                SMCs=fetch_smcs_at_inds(file, list_inds)
                            end

                            #Get all statistics at required times
                            for i in 1:length(times)
                                fork=forks[i]
                                mean_forks[:,i].+=fork

                                mean_zs[:,i].+=fetch_zs(positions[i], fork, N)

                                if !contains(subdir, "No_smcs") && !contains(subdir, "nosmcs") && !contains(subdir, "Ideal")
                                    mean_chipseq[:,i].+=fetch_chipseq(SMCs[i],N)
                                end

                                #fetch_contacts returns all,old,new,inter contacts
                                contacts=fetch_contacts(positions[i], fork, N)
                                for (j,c_array) in enumerate(contacts)
                                    mean_hic[:,:,j,i].+=c_array
                                end

                                av_dist[:,:,i].+=fetch_distance_matrix(positions[i], fork, N)
                            end
                        end
                        num_samp+=1
                end

                #divide by number of samples to get average
                for array in [mean_zs, mean_forks, mean_hic, av_dist]
                    array ./= num_samp
                end

                #save all the statistics to a .h5 file
                out_file_name=replace(dir, data_parent_dir=>out_parent_dir)*".h5"
                h5open(out_file_name, "w") do file
                    @write file num_samp
                    @write file mean_zs
                    @write file mean_forks
                    @write file mean_hic
                    @write file av_dist
                    if !contains(subdir, "No_smcs") && !contains(subdir, "nosmcs") && !contains(subdir, "Ideal")
                        mean_chipseq./=num_samp
                        @write file mean_chipseq
                    end
                end
            end
        end
    end
end

"""
For comparison to Le et al 2013, calculate the normalized Hi-C maps at the same time-points as they have data for
"""
function calc_normalized_hic(times=[0,10,30,45,60,75], N=404,data_parent_dir="Simulation_data/Dynamic/",out_parent_dir="Stats/Dynamic_Hi-C/"; skip_done=true, spacer=1)
    subdirs=readdir(data_parent_dir, join=true)
    subdirs=subdirs[isdir.(subdirs)]

    for subdir in subdirs
        out_dir=replace(subdir, data_parent_dir=>out_parent_dir)
        if !isdir(out_dir)
            mkpath(out_dir)
        end

        simdirs=readdir(subdir, join=true)
        simdirs=simdirs[isdir.(simdirs)]
        for dir in simdirs
            out_file_name=replace(dir, data_parent_dir=>out_parent_dir)*".h5"
            if skip_done && isfile(out_file_name)
                println("Found previous results for $dir, skipping...")
            else
                println("Calculating statistics for $dir")
                num_samp=0
                ter_distinct=true #terminus distinct during dynamic sims
                no_orient=true #for hi-c comparison, don't orient the two strands

                #Initialize output arrays
                mean_forks=zeros(2, length(times))
                mean_hic=zeros(N, N, length(times)) #the 4 is the dimension for all/old/new/inter
                normalized_hic=zeros(N,N, length(times))

                #Loop over files, each containing a single time trajectory
                for filename in readdir(dir,join=true)
                    h5open(filename, "r") do file
                        list_inds=get_indices_at_times(file,times)
                        positions, forks=fetch_pos_fork_at_inds(file, list_inds, no_turn=no_orient, N=N, spacer=spacer, ter_distinct=ter_distinct)

                        #Get all statistics at required times
                        for i in 1:length(times)
                            fork=forks[i]
                            mean_forks[:,i].+=fork

                            #fetch_contacts returns all,old,new,inter contacts
                            contacts=fetch_contacts(positions[i], fork, N)[1]
                            mean_hic[:,:,i].+=contacts
                        end
                    end
                    num_samp+=1
                end

                #divide by number of samples to get average
                for array in [mean_forks, mean_hic]
                    array ./= num_samp
                end

                for time_point in 1:length(times)
                    hic_all=mean_hic[:,:,time_point]
                    hic_all.+=transpose(hic_all)
                    hic_all.*=N/sum(hic_all)
                    normalized_hic[:,:,time_point]=normalized_contact_map(hic_all)
                end

                #save all the statistics to a .h5 file
                out_file_name=replace(dir, data_parent_dir=>out_parent_dir)*".h5"
                h5open(out_file_name, "w") do file
                    @write file num_samp
                    @write file mean_forks
                    @write file mean_hic
                    @write file normalized_hic
                end
            end
        end
    end
end

"""
Calculate time trajectory stats over replication cycle

"""
function calc_time_course_stats(N=404,data_parent_dir="Simulation_data/Dynamic/",out_parent_dir="Stats/Segregation_dynamic/"; skip_done=true, spacer=1)

    subdirs=readdir(data_parent_dir, join=true)
    subdirs=subdirs[isdir.(subdirs)]

    for subdir in subdirs
        out_dir=replace(subdir, data_parent_dir=>out_parent_dir)
        if !isdir(out_dir)
            mkpath(out_dir)
        end

        simdirs=readdir(subdir, join=true)
        simdirs=simdirs[isdir.(simdirs)]
        for dir in simdirs
            out_file_name=replace(dir, data_parent_dir=>out_parent_dir)*"_segregation.h5"
            if skip_done && isfile(out_file_name)
                println("Found previous results for $dir, skipping...")
            else
                println("Calculating segregation statistics for $dir")
                num_samp=zeros(N+1) #number of samples for each stage of replication

                no_orient=!(contains(dir, "no_tether")) #if tethered, oris are distinct

                #Initialize output arrays

                #Track components of the separation between:
                #1. the center of mass of the two replicating regions (strands 1 and 2)
                #2. center of mass of strand 2 and strand 3
                #3. center of mass of strand 1 and strand 3
                #4. the replication forks
                #as a function of R
                z_dist=zeros(4,N+1)
                r_dist=zeros(4,N+1)
                tot_dist=zeros(4,N+1)
                z_dist_std=zeros(4,N+1)
                r_dist_std=zeros(4,N+1)
                tot_dist_std=zeros(4,N+1)

                #segregated fraction of monomers as function of R
                seg_index=zeros(N+1)
                seg_index_std=zeros(N+1)

                #mean long axis position as function of R
                zs=zeros(2*N,N+1)
                zs_std=zeros(2*N,N+1)

                #center of mass long axis position as function of R
                com_zs=zeros(3,N+1)
                com_zs_std=zeros(3,N+1)

                #Loop over files, each containing a single time trajectory
                for filename in readdir(dir,join=true)
                    h5open(filename, "r") do file
                        #Get all statistics at required times
                        for i in keys(file)
                            pos,fork=fetch_pos_fork_at_ind(file,i, no_turn=no_orient, ter_distinct=true, spacer=spacer, N=N)
                            R=fork[1]+N-fork[2]+1

                            num_samp[R]+=1

                            #calculate distances between centers of mass and forks
                            com1,com2,com3=fetch_COMs(pos,fork,N)

                            #center of mass long axis positions
                            com_zs[:,R].+=[com1[3],com2[3],com3[3]]
                            com_zs_std[:,R].+=[com1[3],com2[3],com3[3]].^2

                            #distances between strand centers of mass
                            com_com_12=com1.-com2
                            com_com_23=com3.-com2
                            com_com_13=com1.-com3

                            #distances between forks
                            fork_fork=pos[:,fork[1]].-pos[:,fork[2]]

                            for (index,v) in enumerate([com_com_12, com_com_23, com_com_13, fork_fork])
                                z_dist[index,R]+=abs(v[3])
                                r_dist[index,R]+=norm(v[1:2])
                                tot_dist[index,R]+=norm(v)
                                z_dist_std[index,R]+=abs(v[3])^2
                                r_dist_std[index,R]+=norm(v[1:2])^2
                                tot_dist_std[index,R]+=norm(v)^2
                            end

                            #time lapse of mean long axis positions
                            zs[:,R].+=pos[3,:]
                            zs_std[:,R].+=pos[3,:].^2

                            #segregated fractions
                            si_z=fetch_segregated_fraction(pos,fork,N)

                            @assert si_z<=1 && si_z>=0 "Segregated fraction out of bounds: $si_z !"

                            seg_index[R]+=si_z
                            seg_index_std[R]+=si_z.^2
                        end
                    end
                end

                #divide by number of samples to get average
                for a in [z_dist,r_dist,tot_dist,z_dist_std,r_dist_std,tot_dist_std,zs,zs_std,seg_index,seg_index_std, com_zs,com_zs_std]
                    if ndims(a)>1
                        for i in 1:size(a)[1]
                            a[i,:]./=num_samp
                        end
                    else
                        a./=num_samp
                    end
                end
                
                #calculate standard deviations
                z_dist_std.=sqrt.(round.(z_dist_std.-z_dist.^2, digits=4))
                r_dist_std.=sqrt.(round.(r_dist_std.-r_dist.^2, digits=4))
                tot_dist_std.=sqrt.(round.(tot_dist_std.-tot_dist.^2, digits=4))
                zs_std.=sqrt.(round.(zs_std.-zs.^2, digits=4))
                seg_index_std.=sqrt.(round.(seg_index_std.-seg_index.^2, digits=4))
                com_zs_std.=sqrt.(round.(com_zs_std.-com_zs.^2, digits=4))

                #save all the statistics to a .h5 file
                h5open(out_file_name, "w") do file
                    @write file num_samp
                    @write file r_dist
                    @write file z_dist
                    @write file tot_dist
                    @write file zs
                    @write file com_zs
                    @write file seg_index
                    @write file r_dist_std
                    @write file z_dist_std
                    @write file tot_dist_std
                    @write file zs_std
                    @write file seg_index_std
                    @write file com_zs_std
                end
            end
        end
    end
end

function check_num_samples_dynamic(config_dir="Simulation_data/Dynamic/", out_dir="Stats/")
    folders=vcat(subdirs.(subdirs(config_dir))...)
    num_samples=zeros(Int, length(folders))
    for (index,folder) in enumerate(folders)
        files=readdir(folder, join=true)
        #loop over all files and increment arrays
        for file in files
            h5open(file) do f
                for key in keys(f)
                    i=read_attribute(f[key], "block")+1
                    if i==1
                        num_samples[index]+=1
                    end
                end
            end
        end
    end
    writedlm(out_dir*"num_dynamic_simulations.txt", hcat(replace.(folders, config_dir=>""), num_samples))
end
