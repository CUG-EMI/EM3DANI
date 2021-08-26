#------------------------------------------------------------------------------
#  Define necessary routines for parallel computing.
#------------------------------------------------------------------------------
import Primes.factor

export clear!
export initChannel, initRemoteChannel
export getWorkerIds, assignTasks
export defaultdist, chunk_idxs


#
function initChannel(func::Union{Function,Type}, args::Tuple,kwargs::Array)

    obj = func(args...; kwargs...)
    chan = Channel{typeof(obj)}(1)
    put!(chan, obj)
    return chan

end


# Run the function f with arguments on specified worker.
function initRemoteChannel(func::Union{Function,Type}, pid::Int64, args...; kwargs...)

    return RemoteChannel(()->initChannel(func, args, kwargs), pid)

end


#
function initChannel(data::Any)

    chanRef = Channel{typeof(data)}(1)
    put!(chanRef, data)
    return chanRef

end

function initRemoteChannel(data::Any, pid::Int64=myid())
    rRef = RemoteChannel(()->initChannel(data), pid)
    return rRef

end


# Compute indices array for dividing dims into chunks.
function chunk_idxs(dims, chunks)
    cuts = map(defaultdist, dims, chunks)
    n = length(dims)
    idxs = Array{NTuple{n,UnitRange{Int}}}(undef, chunks...)
    for cidx in CartesianIndices(tuple(chunks...))
        idxs[cidx.I...] = ntuple(i -> (cuts[i][cidx[i]]:cuts[i][cidx[i] + 1] - 1), n)
    end
    return (idxs, cuts)
end


# Get array of start indices for dividing sz into nc chunks.
function defaultdist(sz::Int, nc::Int)
    if sz >= nc
        return ceil.(Int, range(1, stop=sz+1, length=nc+1))
    else
        return [[1:(sz+1);]; zeros(Int, nc-sz)]
    end
end

# Decide how to divide each dimension, returns size of chunks array
function defaultdist(dims, pids)
    dims = [dims...]
    chunks = ones(Int, length(dims))
    np = length(pids)
    f = sort!(collect(keys(factor(np))), rev=true)
    k = 1
    while np > 1
        # repeatedly allocate largest factor to largest dim
        if np % f[k] != 0
            k += 1
            if k > length(f)
                break
            end
        end
        fac = f[k]
        (d, dno) = findmax(dims)
        # resolve ties to highest dim
        dno = findlast(isequal(d), dims)
        if dims[dno] >= fac
            dims[dno] = div(dims[dno], fac)
            chunks[dno] *= fac
        end
        np = div(np, fac)
    end
    return chunks
end


function initRemoteChannel(data::Vector, pids::Vector{Int64})

    np = length(pids)
    pRef = Array{RemoteChannel}(undef, np)
    nd   = length(data)
    ival = defaultdist(nd, np)
    for (idx,ip) in enumerate(pids)
        pRef[idx] = initRemoteChannel(data[ival[idx]:ival[idx+1]-1], ip)
    end
    return pRef

end


function initRemoteChannel(A::AbstractArray;
    procs = workers()[1:min(nworkers(), maximum(size(A)))],
    dist = defaultdist(size(A), procs))

    #
    np = prod(dist)
    procs_used = procs[1:np]
    idxs, _ = chunk_idxs([size(A)...], dist)
    pRef = Array{RemoteChannel}(undef, np)
    for (idx,ip) in enumerate(procs_used)
        pRef[idx] = initRemoteChannel(A[idxs[idx]...], ip)
    end

    return pRef

end


#
function clear!(r::RemoteChannel)

    p = take!(r)
    p = clear!(p)
    put!(r, p)

    return r
end


function clear!(f::Future)

    p = take(f)
    p = clear!(p)
    put!(f, p)

    return f
end


function clear!(pf::Array{RemoteChannel})

    @sync begin
        for p=workers()
            @async begin
                for i=1:length(pf)
                    if p == pf[i].where
                        pf[i] = remotecall_fetch(clear!, p, pf[i])
                    end
                end
            end
        end
    end

end


function clear!(pf::Array{Future})

    @sync begin
        for p=workers()
            @async begin
                for i=1:length(pf)
                    if p == pf[i].where
                        pf[i] = remotecall(clear!, p, pf[i])
                    end
                end
            end
        end
    end

end


#
function getWorkerIds(A::Array{RemoteChannel})
    Ids = []
    for k=1:length(A)
        push!(Ids,A[k].where)
    end
    return unique(Ids)
end


# Distribute tasks over processors as evenly as possible.
function assignTasks(nTask::Int, nProc::Int, method::Int=1)

    nProc > nTask  &&  error("There are more processors than tasks!")

    idx = Array{Any}(undef, nProc)

    if method == 1
        k = fld(nTask, nProc)
        r = rem(nTask, nProc)

        for i=1:r
            id1 = (i-1)*(k+1) + 1
            id2 = id1 + k
            idx[i] = id1:id2
        end

        for i=r+1:nProc
            id1 = (i-1)*k + 1 + r
            id2 = id1 + k-1
            idx[i] = id1:id2
        end

    elseif method == 2
        for i = 1:nProc
            idx[i] = i:nProc:nTask
        end

    end

    return idx
end
