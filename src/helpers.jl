"""
File contains functions for parsing, loading simulation files, and converting between
replicated length, time, etc.

Loop-extruder simulations save .h5 files that store the following for each time-point i:
i/pos   ->  coordinates of monomers, array [3, N]
i/block ->  simulation step
i/Fork  ->  current coordinates of replication fork, array [2]
i/SMCs  ->  positions of loop-extruder links, array [2,M] for non-topological, 
            for topological the number of entries depends on the 
            number of strands in each loop-extruder leg
This file provides functions for fetching data for a given set of simulations.
"""

using HDF5

periodic_ind(i,N)= i>1 ? Int((i-1)%N+1) : periodic_ind(i+N,N)

replicated(fork,N) = vcat(fork[2]+1:N, 1:fork[1]-1)

unreplicated(fork)=fork[1]:fork[2]

function parse_after(dir_name, what="/R_")
    from=findfirst(what, dir_name)[end]+1
    to=from+findfirst("_", dir_name[from:end])[end]-2
    R=parse(Int,dir_name[from:to])
    return R
end

function fork_to_R(fork,N)
    return sort(fork)[1]+N-sort(fork)[2]
end

function R_to_time(n_kb; r_replication=53.9) #r_replication is speed of replication in kb/min
    return n_kb/r_replication
end

function fork_to_time(fork; N=404, mon_kb=10, r_replication=53.9) #r_replication is speed of replication in kb/min
    n_kb=fork_to_R(fork,N)*mon_kb
    return n_kb/r_replication
end

time_to_block(time, blocks_per_minute=0.559)=round(Int,time/blocks_per_minute)
block_to_time(block, blocks_per_minute=0.559)=block*blocks_per_minute

time_to_height(t,r_0=25.6*88/1000, growth_exponent=0.0055)=r_0*exp(growth_exponent*t)

function fork_to_height(fork; N=404, mon_kb=10, r_0=25.6*88/1000, growth_exponent=0.0055)
    t=fork_to_time(fork, N=N, mon_kb=mon_kb)
    return r_0*exp(growth_exponent*t)
end

function R_to_height(n; mon_kb=10, r_replication=53.9, r_0=25.6*88/1000, growth_exponent=0.0055)
    t=n*mon_kb/r_replication
    return r_0*exp(growth_exponent*t)
end


#Shift map so that the middle of the map is at the middle of the matrix
function shifted_map(M,mid=1)
    if isnan(mid)
        return M
    end
    N=size(M)[1]
    if N%2==0
            first=Int(mid-N/2)
            last=Int(mid+N/2-1)
    else
            first=Int(mid-((N-1)/2))
            last=Int(mid+(N-1)/2)
    end
    inds=periodic_ind.(first:last,N)
    return M[inds,inds]
end

function shifted_map_replicated(M,mid=1)
    if isnan(mid)
        return M
    end
    N=Int(size(M)[1]/2)
    M_shifted=copy(M)
    for inds in [1:N, N+1:2*N]
        for sec_inds in [1:N, N+1:2*N]
            M_shifted[inds, sec_inds].=shifted_map(M[inds,sec_inds],mid)
        end
    end
    return M_shifted
end

function shifted_vector(v,mid=1)
    N=size(v)[1]
    if N%2==0
            first=Int(mid-N/2)
            last=Int(mid+N/2-1)
    else
            first=Int(mid-((N-1)/2))
            last=Int(mid+(N-1)/2)
    end
    inds=periodic_ind.(first:last,N)
    return v[inds]
end

function shifted_vector_replicated(v,mid=1)
    N=Int(size(v)[1]/2)
    if N%2==0
            first=Int(mid-N/2)
            last=Int(mid+N/2-1)
    else
            first=Int(mid-((N-1)/2))
            last=Int(mid+(N-1)/2)
    end
    inds=vcat(periodic_ind.(first:last,N),periodic_ind.(first:last,N).+N)
    return v[inds]
end

#return a list of all the subdirectories
function subdirs(dir::String)
    sdirs=readdir(dir, join=true)
    return sdirs[isdir.(sdirs)]
end

function fetch_time(n, file)
    read_attribute(file["$n"], "block")
    return block_to_time(block)
end

"""
Switch monomer indices if arm 1 is further from pole 1
"""
function orient_arms_chromosome!(array,fork,N=404)
    @assert size(array)[2]==2*N
    ter=round(Int,N/2)
    #check if need to relabel arms
    if mean(array[3,1:ter])>mean(array[3,ter:N])
        #arm2 closer to pole 1; reverse the indices
        array[:,1:N].=reverse(array[:,1:N], dims=2)
        array[:,N+1:2*N].=reverse(array[:,N+1:2*N], dims=2)
        fork=reverse(N .- fork)
    end
end

"""
Switch strand labels if center of mass of strand 2 is closer to a pole

Also flip so that pole 1 is at -L/2
"""
function com_turn_chromosome!(array, fork, N=404, ter_distinct=true)
    @assert size(array)[2]==2*N
    replicated=vcat(1:fork[1], fork[2]:N)
    replicated_2=replicated.+N
    unreplicated=vcat(floor(Int,N/2):-1:fork[1]+1, fork[2]-1:-1:floor(Int, N/2)+1) #ter to fork1, fork2 to ter
    ter=round(Int,N/2)

    if !ter_distinct && length(unreplicated)==length(replicated)
        com1=mean(array[3,replicated])
        com2=mean(array[3,replicated_2])
        com3=mean(array[3,unreplicated])
        coms=[com1,com2,com3]
        order=sortperm(coms)
        sort!(coms)
        if abs(coms[1])<abs(coms[3])
            reverse!(order)
        end

        strand1, strand2, strand3 = [replicated, replicated_2, unreplicated][order]
        array[:, replicated].=array[:, strand1]
        array[:, unreplicated].=array[:, strand2]
        array[:, replicated_2].=array[:, strand3]

    else
        #one strand different length than others
        com1=mean(array[3,1:N])
        com2=mean(array[3,vcat(replicated_2,unreplicated)])        
        #Strand 1 should have a center of mass closer to mid-cell, ie 0
        if abs(com1)>abs(com2)
            swapped=copy(array)
            #switch labels old and new
            array[:,replicated].=swapped[:,replicated.+N]
            array[:,replicated.+N].=swapped[:,replicated]
            com1=com2
        end
    end
    #Strand 1 should always be near negative pole
    if mean(array[3,replicated])>0
        array[3,:].*=-1
    end

    #check if need to relabel arms
    if mean(array[3,1:ter])>mean(array[3,ter:N])
        #arm2 closer to pole 1; reverse the indices
        array[:,1:N].=reverse(array[:,1:N], dims=2)
        array[:,N+1:2*N].=reverse(array[:,N+1:2*N], dims=2)
        fork=reverse(N .- fork)
    end
end

"""
This function calculates the time for each index based on the block number,
returns the indices for the time-points closest to those in array
'list_times'
"""
function get_indices_at_times(file, list_times)
    ks=keys(file)
    blocks=map(k-> read_attribute(file["$k"], "block"), ks)
    required_blocks=time_to_block.(list_times)
    list_inds=[]
    for required_block in required_blocks
        i=findfirst(blocks .== required_block)
        push!(list_inds, ks[i])
    end
    return list_inds
end

function get_simulation_step(file, index)
    t=read_attribute(file["$index"], "block")
    return t
end

"""
Get the array for monomer positions
"""
function fetch_pos_at_inds(file, inds, R; no_turn=true, N=404, spacer=1, ter_distinct=true)
    data=[]
    for i in inds
        pos=read(file["$i/pos"])[:,1:spacer:2*spacer*N]
        fork=round.(Int,[R/spacer/2,N-R/spacer/2])
        if no_turn
            orient_arms_chromosome!(pos,fork,N)
        else
            com_turn_chromosome!(pos,fork,N, ter_distinct)
        end
        push!(data,pos)
    end
    return data
end

"""
Fetch position array for a single index
"""
function fetch_pos_at_ind(file, ind, R; no_turn=true, N=404, spacer=1, ter_distinct=true)
    pos=read(file["$ind/pos"])[:,1:spacer:2*spacer*N]
    fork=round.(Int,[R/spacer/2,N-R/spacer/2])

    if no_turn
        orient_arms_chromosome!(pos,fork,N)
    else
        com_turn_chromosome!(pos,fork,N, ter_distinct)
    end
    return pos
end

"""
Fetch position array for a single index
"""
function fetch_raw_pos_at_ind(file, ind; N=404, spacer=1)
    pos=read(file["$ind/pos"])[:,1:spacer:2*spacer*N]
    return pos
end

"""
Get the array for monomer positions
"""
function fetch_pos_fork_at_inds(file, inds; no_turn=true, N=404, spacer=1, ter_distinct=true)
    data=[]
    forks=[]
    for i in inds
        pos=read(file["$i/pos"])[:,1:spacer:2*spacer*N]
        fork=fetch_fork_at_ind(file,i)
        if no_turn
            orient_arms_chromosome!(pos,fork,N)
        else
            com_turn_chromosome!(pos,fork,N,ter_distinct)
        end
        push!(data,pos)
        push!(forks,fork)
    end
    return data, forks
end

"""
Get the array for monomer positions without relabelling them; good for animations
"""
function fetch_raw_pos_fork_at_inds(file, inds; N=404, spacer=1)
    data=[]
    forks=[]
    ter=floor(Int, N/2)
    for i in inds
        pos=read(file["$i/pos"])[:,1:spacer:2*spacer*N]
        fork=fetch_fork_at_ind(file,i)
        push!(data,pos)
        push!(forks,fork)
    end
    return data, forks
end

"""
Fetch position array and fork for a single index
"""
function fetch_pos_fork_at_ind(file, ind; no_turn=true, N=404, spacer=1, ter_distinct=true)
    pos=read(file["$ind/pos"])[:,1:spacer:2*spacer*N]
    fork=fetch_fork_at_ind(file,ind)
    if no_turn
        orient_arms_chromosome!(pos,fork,N)
    else
        com_turn_chromosome!(pos,fork,N,ter_distinct)
    end
    return pos, fork
end

"""
For a set of indices, fetch fork positions/spacer
from .h5 file

These can then be used as indices for replicated and unreplicated
sections of the chromosome
"""
function fetch_fork_at_inds(file, inds; spacer=1)
    data=[]
    for i in inds
        push!(data,floor.(Int,read(file["$i/Fork"])./spacer).+1)
    end
    return data
end

function fetch_fork_at_ind(file, ind; spacer=1)
    return floor.(Int,read(file["$ind/Fork"])./spacer.+1)
end

"""
Given .h5 file, fetch the positions of SMC links along chromosome
"""
function fetch_smcs_at_inds(file, inds; spacer=1)
    data=[]
    for i in inds
        push!(data,floor.(Int,read(file["$i/SMCs"])./spacer).+1)
    end
    return data
end

function fetch_smcs_at_ind(file, ind; spacer=1)
    return floor.(Int,read(file["$ind/SMCs"])./spacer).+1
end

"""
Given a contact map, return the normalized contact map.

Normalization as for Le et al. 2013: each row and column sum up to 1
"""
function normalized_contact_map(hi_c_map, num_steps=200)
    normalized=copy(hi_c_map)
    N=size(hi_c_map)[1]

    for k in 1:num_steps
        #one normalisation step
        total_reads=sum(normalized)
        totals=zeros(N)
        for i in 1:N
            totals[i]=sum(normalized[:,i])
        end
        for i in 1:N
            for j in 1:N
                normalized[j,i]*=total_reads/(totals[i]*totals[j])/N
            end
        end
    end
    normalized.*=N/sum(normalized)

    return normalized
end

"""
Function that returns matrix with A in lower half
and B in upper half. For plotting two matrices at once.
"""
function half_half(A,B)
    N=size(A)[1]
    halves=zeros(N,N)
    for i in 1:N
        for j in i+1:N
            halves[i,j]=A[i,j]
            halves[j,i]=B[i,j]
        end
    end
    return halves
end
