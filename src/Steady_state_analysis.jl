
include("Statistics_functions.jl")
using DelimitedFiles

"""
Given a directory of configurations performed at a fixed replication stage, calculate...

- mean long axis positions
- segregated fractions
- fork and center of mass separations

averaged over configurations from simulation step sample_from onwards
"""
function calc_fast_steady_state_stats(config_dir,out_dir="Stats/Equilibrium_stats/"; skip_done=true, sample_from=0, spacer=1, N=404)
    for subdir in subdirs(config_dir)
        for dir in subdirs(subdir)
            out_path=replace(dir, config_dir=>out_dir)
            if !isfile(out_path*"steady_state_stats.h5") || !skip_done
                #read extent of replication from directory name
                R=parse_after(dir)
                no_orient=!(contains(dir, "no_tether")) #if tethered, dont relabel oris
                ter_distinct=!(contains(dir, "No_smcs")) #if smcs, ter distinct from ori distinct even if R=N/2

                println("Calculating average z-positions and segregated fractions for $dir with R=$R")
                fork=round.(Int,[R/spacer/2,N-R/spacer/2])

                #arrays for stats
                num_samp=0
                mean_zs=zeros(2*N)
                zs_std=zeros(2*N)

                #segragation stats
                z_dist=zeros(4)
                r_dist=zeros(4)
                tot_dist=zeros(4)
                z_dist_std=zeros(4)
                r_dist_std=zeros(4)
                tot_dist_std=zeros(4)

                seg_index=0
                seg_index_std=0

                com_zs=zeros(3)
                com_zs_std=zeros(3)

                for f in readdir(dir, join=true)
                    h5open(f) do file
                        for i in keys(file)
                            if read_attribute(file[i], "block")>sample_from
                                pos=fetch_pos_at_ind(file, i, R, no_turn=no_orient,N=N,spacer=spacer, ter_distinct=ter_distinct)
                                z=fetch_zs(pos,fork, N)
                                mean_zs.+=z
                                zs_std.+=z.^2

                                seg_index+=fetch_segregated_fraction(pos,fork,N)
                                seg_index_std+=fetch_segregated_fraction(pos,fork,N)^2

                                #segregation stats
                                com1,com2,com3=fetch_COMs(pos,fork,N)

                                com_zs.+=[com1[3],com2[3],com3[3]]
                                com_zs_std.+=[com1[3],com2[3],com3[3]].^2
                                #distances
                                com_com_12=com1.-com2
                                com_com_23=com3.-com2
                                com_com_13=com1.-com3

                                fork_fork=pos[:,fork[1]].-pos[:,fork[2]]

                                for (index,v) in enumerate([com_com_12, com_com_23, com_com_13, fork_fork])
                                    z_dist[index]+=abs(v[3])
                                    r_dist[index]+=norm(v[1:2])
                                    tot_dist[index]+=norm(v)
                                    z_dist_std[index]+=abs(v[3])^2
                                    r_dist_std[index]+=norm(v[1:2])^2
                                    tot_dist_std[index]+=norm(v)^2
                                end

                                num_samp+=1
                            end
                        end
                    end
                end

                if num_samp>0
                    #divide by number of samples to get average
                    for array in [mean_zs, com_zs, com_zs_std,
                                    z_dist, r_dist, tot_dist,z_dist_std,
                                    r_dist_std,tot_dist_std, zs_std]
                        array ./= num_samp
                    end
                    seg_index/=num_samp
                    seg_index_std/=num_samp

                    #calculate stds
                    seg_index_std=sqrt(seg_index_std-seg_index^2)
                    zs_std.=sqrt.(zs_std.-mean_zs.^2)
                    z_dist_std.=sqrt.(z_dist_std.-z_dist.^2)
                    r_dist_std.=sqrt.(r_dist_std.-r_dist.^2)
                    tot_dist_std.=sqrt.(tot_dist_std.-tot_dist.^2)
                    com_zs_std.=sqrt.(com_zs_std.-com_zs.^2)

                    #save to file
                    if !ispath(out_path)
                        mkpath(out_path)
                    end
                    h5open(out_path*"/steady_state_stats.h5", "w") do file
                        @write file num_samp
                        @write file mean_zs
                        @write file zs_std
                        @write file fork
                        @write file com_zs
                        @write file com_zs_std
                        @write file r_dist
                        @write file z_dist
                        @write file tot_dist
                        @write file r_dist_std
                        @write file z_dist_std
                        @write file tot_dist_std
                        @write file seg_index
                        @write file seg_index_std
                    end
                end
            else
                println("Found results in $out_path, skipping calculations")
            end
        end
    end
end

"""
Given a directory of configurations performed at a fixed replication stage, calculate...

- contact probability maps
- distance maps

averaged over configurations from simulation step sample_from onwards
"""
function calc_slow_steady_state_stats(config_dir,out_dir="Stats/Equilibrium_stats/"; skip_done=true, sample_from=0)
    for subdir in subdirs(config_dir)
        for dir in subdirs(subdir)
            out_path=replace(dir, config_dir=>out_dir)
            if !isfile(out_path*"/slow_steady_state_stats.h5") || !skip_done
                no_orient=!(contains(dir, "no_tether")) #if tethered, dont relabel oris
                ter_distinct=!(contains(dir, "No_smcs")) #if smcs, ter distinct from ori distinct even if R=N/2
                
                #read extent of replication from directory name
                R=parse_after(dir)
                N=parse_after(dir, "_N_")
                spacer=round(Int,N/405)
                N=round(Int, N/spacer)

                println("Calculating distance and Hi-C maps for dir $dir with R=$R")
                fork=round.(Int,[R/spacer/2,N-R/spacer/2])

                #arrays for stats
                num_samp=0
                mean_hic=zeros(N, N, 4) #the 4 is the dimension for all/old/new/inter
                av_dist=zeros(2*N,2*N)

                for f in readdir(dir, join=true)
                    h5open(f) do file
                        for i in keys(file)
                            if read_attribute(file[i], "block")>sample_from
                                pos=fetch_pos_at_ind(file, i, R, no_turn=no_orient,N=N,spacer=spacer, ter_distinct=ter_distinct)

                                #fetch_contacts returns all,old,new,inter contacts
                                contacts=fetch_contacts(pos, fork, N)
                                for (j,c_array) in enumerate(contacts)
                                    mean_hic[:,:,j].+=c_array
                                end

                                av_dist.+=fetch_distance_matrix(pos, fork, N)

                                num_samp+=1
                            end
                        end
                    end
                end

                if num_samp>0
                    #divide by number of samples to get average
                    for array in [av_dist,mean_hic]
                        array ./= num_samp
                    end

                    #save to file
                    if !ispath(out_path)
                        mkpath(out_path)
                    end
                    h5open(out_path*"/slow_steady_state_stats.h5", "w") do file
                        @write file num_samp
                        @write file fork
                        @write file mean_hic
                        @write file av_dist
                    end

                    #center of mass array for calculating statistics
                    coms=zeros(3, num_samp)
                    counter=1
                    for f in readdir(dir, join=true)
                        h5open(f) do file
                            for i in keys(file)
                                if read_attribute(file[i], "block")>sample_from
                                    pos=fetch_pos_at_ind(file, i, R, no_turn=no_orient,N=N,spacer=spacer, ter_distinct=ter_distinct)
                                    coms[:,counter].+=fetch_z_COMs(pos,fork, N)
                                    counter+=1
                                end
                            end
                        end
                    end

                    com_means=mean(coms, dims=2)
                    com_covariances=cov(coms, dims=2)

                    writedlm(out_path*"com_covariances.txt", hcat(com_means,com_covariances))
                end
            else
                println("Found results in $out_path, skipping calculations")
            end
        end
    end
end

"""
Given a directory of configurations performed at a fixed replication stage, calculate...

- mean long axis positions
- center of mass positions
- segregated fractions

as a function of simulation step
"""
function check_convergence(config_dir="Simulation_data/Steady_state/", out_dir="Stats/Convergence/"; max_it=2500, N=404, spacer=1)
    folders=vcat(subdirs.(subdirs(config_dir))...)

    if !isdir(out_dir*"Steady_state_No_smcs")
        mkpath(out_dir)
        for subdir in ["Steady_state_No_smcs", "Steady_state_No_smcs_no_tether", "Steady_state_Nontopological_smcs_no_tether", "Steady_state_Nontopological_smcs", "Steady_state_Topological_smcs_no_tether", "Steady_state_Topological_smcs"]
            mkpath(out_dir*subdir)
        end
    end

    for folder in folders
        files=readdir(folder, join=true)
        no_orient=!(contains(folder, "no_tether")) #if tethered, dont relabel oris
        ter_distinct=!(contains(folder, "No_smcs")) #if smcs, ter distinct from ori distinct even if R=N/2

        zs=zeros(2*N,max_it)
        coms=zeros(3, max_it)
        seg_inds=zeros(max_it)
        counts=zeros(max_it)
        
        #loop over all files and increment arrays
        for file in files
            R=parse_after(file)
            fork=Int.([R/2/spacer, N-R/2/spacer])
            h5open(file) do f
                for key in keys(f)
                    i=read_attribute(f[key], "block")+1
                    if i>max_it
                        println(file)
                        println("got $i")
                    else
                        pos=fetch_pos_at_ind(f,key,R, spacer=spacer, N=N, no_turn=no_orient, ter_distinct=ter_distinct)
                        zs[:,i].+=fetch_zs(pos,fork,N)
                        coms[:,i].+=fetch_z_COMs(pos,fork,N)
                        seg_inds[i]+=fetch_segregated_fraction(pos,fork,N)
                        counts[i]+=1
                    end
                end
            end
        end

        #divide by counts to get averages
        for r in eachrow(zs)
            r./=counts
        end
        for r in eachrow(coms)
            r./=counts
        end
        seg_inds./=counts
        
        #save the data
        out_path=replace(folder, config_dir=>out_dir)
        writedlm(out_path*"z_over_sims.txt", zs)
        writedlm(out_path*"com_over_sims.txt", coms)
        writedlm(out_path*"seg_ind_over_sims.txt", seg_inds)
    end
end