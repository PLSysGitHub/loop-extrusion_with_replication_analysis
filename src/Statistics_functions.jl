include("helpers.jl")

#helper functions for fetching time point data from time trajectory simulations
using LinearAlgebra, StatsBase, Distances

function fetch_zs(pos, fork, N=405)
    zs=pos[3,:]
    #unreplicated; just use coordinate of first
    zs[N+fork[1]+1:N+fork[2]-1].=pos[3,fork[1]+1:fork[2]-1]
    return zs
end

function fetch_z_COMs(pos, fork, N=404)
    zs=pos[3,:]
    #unreplicated; just use coordinate of first
    return mean(zs[replicated(fork,N)]), mean(zs[replicated(fork,N).+N]), mean(zs[unreplicated(fork)])
end

function fetch_COMs(pos, fork, N=404)
    #unreplicated; just use coordinate of first
    return mean(pos[:,replicated(fork,N)], dims=2), mean(pos[:,replicated(fork,N).+N],dims=2), mean(pos[:,unreplicated(fork)], dims=2)
end

function fetch_contacts(pos,fork, N=404, contact_r=2)
    contacts_old=zeros(N,N)
    contacts_new=zeros(N,N)
    contacts_inter=zeros(N,N)
    repl=union(N+1:max(N+1,fork[1]+N), min(2*N,fork[2]+N):2*N)
    #old-old
    for i in 1:N
        for j in i+2:N
            if norm(pos[:,i].-pos[:,j])<contact_r
                contacts_old[j,i]+=1
            end
        end
    end
    #new-new
    for i in repl
        for j in repl[repl.>i+1]
            if norm(pos[:,i].-pos[:,j])<contact_r
                contacts_new[j-N,i-N]+=1
            end
        end
    end
    #inter
    for i in 1:N
        for j in repl
            if abs(i-periodic_ind(j,N))>1
                if norm(pos[:,i].-pos[:,j])<contact_r
                    contacts_inter[j-N,i]+=1
                end
            end
        end
    end
    contacts=contacts_old.+contacts_new.+contacts_inter
    return contacts, contacts_old, contacts_new, contacts_inter
end

function fetch_chipseq(smcs, N=404)
    inds=reshape(smcs, :,1)
    inds[inds.>N].-=N
    chipseq=zeros(N)
    chipseq[inds].+=1
    return chipseq
end

function fetch_distance_matrix(pos, fork, N=404)
    not_replicated=N+fork[1]+1:N+fork[2]-1

    d_matrix =pairwise(Euclidean(), pos, pos, dims=2)
    d_matrix[not_replicated,1:N].=d_matrix[not_replicated.-N,1:N]
    d_matrix[1:N,not_replicated].=d_matrix[1:N,not_replicated.-N]
    d_matrix[not_replicated,N+1:2*N].=d_matrix[not_replicated.-N,1:N]
    d_matrix[N+1:2*N,not_replicated].=d_matrix[1:N,not_replicated.-N]
    return d_matrix
end

function fetch_segregated_fraction(pos, fork, N=404)
    R=fork[1]+N-fork[2]
    zs=pos[3,:]
    below=length(findall(zs[replicated(fork,N)].<0))
    above=length(findall(zs[replicated(fork,N).+N].>0))
    return abs(above+below-R)/R
end

function fetch_mean_squared_dist_array(pos,N=404, linear=true)
    ds=pairwise(Euclidean(), pos, dims=2) #efficient distance matrix calculation

    if linear
        mean_sq=zeros(N-1)
        counts=zeros(N-1)
        for i in 1:N
            for j in i+1:N
                s=j-i
                mean_sq[s]+=ds[j,i]^2
                counts[s]+=1
            end
        end
        mean_sq./=counts
    else
        mean_sq=zeros(floor(Int,N/2))
        counts=zeros(floor(Int,N/2))
        for i in 1:N
            for j in i+1:N
                s=min(j-i, N-(j-i))
                mean_sq[s]+=ds[j,i]^2
                counts[s]+=1
            end
        end
        mean_sq./=counts
    end
    return mean_sq
end