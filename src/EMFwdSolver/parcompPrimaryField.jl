export parPrimaryElectricField, compPrimaryEfieldFreqTx
#-------------------------------------------------------------------------------
"""
    compPrimaryEfieldFreqTx(emMesh, emData, refModel, freqTxIdx)

calculates primary field at specific frequency and transmitter.

"""
function compPrimaryEfieldFreqTx(emMesh::EMTensorMesh, emData::CSEMData,
                                 refModel::RefModel, freqTxIdx::Array)

    #
    dipLen = emData.dpLen

    # grid location
    xNode = cumsum([0; emMesh.xLen]) .- emMesh.origin[1]
    yNode = cumsum([0; emMesh.yLen]) .- emMesh.origin[2]
    zNode = cumsum([0; emMesh.zLen]) .- emMesh.origin[3]

    # get sampling locations of EM fields
    # first get cell center coordiante
    xCen = xNode[1:end-1] + diff(xNode) / 2
    yCen = yNode[1:end-1] + diff(yNode) / 2
    zCen = zNode[1:end-1] + diff(zNode) / 2

    # get background model
    sig1D   = refModel.sig1D
    depth1D = refModel.depth1D

    # exclude air layer and seawater layer
    nAir = length(emMesh.airLayer)
    nSea = length(emMesh.seaLayer)

    #
    zlayer = nAir + nSea
    zNode  = zNode[zlayer+1:end]
    zCen   = zCen[zlayer+1:end]

    # edges and faces
    nx, ny, nz = emMesh.gridSize
    nEx = nx * (ny+1) * zlayer
    nEy = (nx+1) * ny * zlayer
    nEz = (nx+1) * (ny+1) * zlayer
    pEx = zeros(ComplexF64, nEx)
    pEy = zeros(ComplexF64, nEy)
    pEz = zeros(ComplexF64, nEz)

    # edge center
    tmpx, tmpy, tmpz = ndgrid(xCen, yNode, zNode)
    exLoc = [tmpx[:] tmpy[:] tmpz[:] ]
    tmpx, tmpy, tmpz = ndgrid(xNode, yCen, zNode)
    eyLoc = [tmpx[:] tmpy[:] tmpz[:] ]
    tmpx, tmpy, tmpz = ndgrid(xNode, yNode, zCen)
    ezLoc = [tmpx[:] tmpy[:] tmpz[:] ]


    # calculate primary field
    epField = zeros(ComplexF64, 0)   # will have size (nFwd*nE, 1)
    nFwd = size(freqTxIdx, 1)
    for i=1:nFwd
        txIdx, freqIdx = freqTxIdx[i]
        println("Calculating primary field at freqID=$(freqIdx), txID=$(txIdx) ...")

        freq = emData.freqArray[freqIdx]
        freq = collect(freq)
        rxtmpz = zeros(1, 3)
        txLoc  = zeros(1, 5)
        txLoc[1, :] = emData.txLoc[txIdx, :]

        if dipLen < 10.0
            emFlag = 1
            @printf("%-12s", "     - Ex:")
            t = @elapsed ex1d = comp1DCSEM(txLoc, exLoc, freq, sig1D, depth1D, emFlag)
            @printf("%8.2f %s\n", t, "seconds.")
            ex1d = vcat(pEx, vec(ex1d))

            emFlag = 2
            @printf("%-12s", "     - Ey:")
            t = @elapsed ey1d = comp1DCSEM(txLoc, eyLoc, freq, sig1D, depth1D, emFlag)
            @printf("%8.2f %s\n", t, "seconds.")
            ey1d = vcat(pEy, vec(ey1d))

            emFlag = 3
            @printf("%-12s", "     - Ez:")
            ez1d = zeros(ComplexF64, 0)
            t = @elapsed begin
                @inbounds for k = 1:size(ezLoc, 1)
                    rxtmpz[1, :] = ezLoc[k, :]
                    tmp = comp1DCSEM(txLoc, rxtmpz, freq, sig1D, depth1D, emFlag)
                    append!(ez1d, collect(vec(tmp)))
                end
            end
            @printf("%8.2f %s\n", t, "seconds.")
            ez1d = vcat(pEz, vec(ez1d))

        else # finite length
            println("Finite length source is used!")
            nIntPts = 15
            compFlag = 1
            @printf("%-12s", "     - Ex:")
            t = @elapsed ex1d = dipLen * comp1DCSEM(txLoc, exLoc, dipLen, nIntPts, freq,
                                             sig1D, depth1D, compFlag)
            @printf("%8.2f %s\n", t, "seconds.")
            ex1d = vcat(pEx, vec(ex1d))

            compFlag = 2
            @printf("%-12s", "     - Ey:")
            t = @elapsed ey1d = dipLen * comp1DCSEM(txLoc, eyLoc, dipLen, nIntPts, freq,
                                             sig1D, depth1D, compFlag)
            @printf("%8.2f %s\n", t, "seconds.")
            ey1d = vcat(pEy, vec(ey1d))

            compFlag = 3
            @printf("%-12s", "     - Ez:")
            ez1d = zeros(ComplexF64, 0)
            t = @elapsed begin
                @inbounds for k = 1:size(ezLoc, 1)
                    rxtmpz[1, :] = ezLoc[k, :]
                    tmp = dipLen * comp1DCSEM(txLoc, rxtmpz, dipLen, nIntPts, freq,
                                              sig1D, depth1D, compFlag)
                    append!(ez1d, collect(tmp))
                end
            end
            @printf("%8.2f %s\n", t, "seconds.")
            ez1d = vcat(pEz, vec(ez1d))

        end

        append!(epField , vcat(ex1d, ey1d, ez1d))

    end

    return epField

end


#
# This routine is particularly for MT, while the above one is for CSEM.
#
function compPrimaryEfieldFreqTx(emMesh::EMTensorMesh, emData::MTData,
                                 refModel::RefModel, freqTxIdx::Array)

    #
    sig1D   = refModel.sig1D
    depth1D = refModel.depth1D
    nx, ny, nz = emMesh.gridSize

    zNode = cumsum([0; emMesh.zLen]) .- emMesh.origin[3]

    sigLayer = ones(nz) * sig1D[end]
    for j=1:length(depth1D)-1
        idx1 = findfirst(isequal(depth1D[j]), zNode)
        idx2 = findfirst(isequal(depth1D[j+1]), zNode)
        sigLayer[idx1:idx2-1] .= sig1D[j]
    end

    nEx = nx * (ny+1) * (nz+1)
    nEy = (nx+1) * ny * (nz+1)
    nEz = (nx+1) * (ny+1) * nz
    nE = nEx + nEy + nEz
    pEx = zeros(ComplexF64, nEx)
    pEy = zeros(ComplexF64, nEy)
    pEz = zeros(ComplexF64, nEz)


    # calculate primary field
    epField = zeros(ComplexF64, 0)   # will have size (nFwd*nE, 1)
    nFwd = size(freqTxIdx, 1)
    for i=1:nFwd
        txIdx, freqIdx = freqTxIdx[i]
        println("Calculating primary field at freqID=$(freqIdx), txID=$(txIdx) ...")

        freq = emData.freqArray[freqIdx]
        epTmp = mt1DAnalyticFields(freq, sigLayer, zNode)

        if txIdx == 1
            Ex = zeros(ComplexF64, nx, ny+1, nz+1)
            for iz = 1:nz+1
                Ex[:, :, iz] .= epTmp[iz]
            end
            e1d = vcat(vec(Ex), pEy, pEz)

        elseif txIdx == 2
            Ey = zeros(ComplexF64, nx+1, ny, nz+1)
            for iz = 1:nz+1
                Ey[:, :, iz] .= epTmp[iz]
            end
            e1d = vcat(pEx, vec(Ey), pEz)

        end

        append!(epField , e1d)

    end

    return epField

end


#
"""
    parPrimaryElectricField(emMesh, emData, refModel, pids)

calculate primary field parallelly.

"""
function parPrimaryElectricField(emMesh::EMTensorMesh, emData::EMData,
                                 refModel::RefModel, pids::Vector{Int})

    #
    if typeof(emData) <: CSEMData
        probType = "csem"
        ns = size(emData.txLoc, 1)
    elseif typeof(emData) <: MTData
        probType = "mt"
        ns = 2
    end

    nFreq = length(emData.freqArray)
    nFwd  = nFreq * ns

    np = length(pids)
    if np > nFwd
        @warn("There are more processors than tasks.")
        rmprocs(pids[nFwd+1:end])
        pids = pids[1:nFwd]
    else
        k = fld(nFwd, np)
        printstyled("there are $(k) task(s) assigned to each processor.\n", color=:cyan)
    end
    np = length(pids)

    #
    fwdIdx = assignTasks(nFwd, np)

    # get the (txIdx, freqIdx) pairs
    freqTxIdx = Array{Tuple}(undef, nFwd)
    m=0
    for i=1:nFreq
        for j=1:ns
            m += 1
            freqTxIdx[m] = (j, i)
        end
    end

    #
    epMat = Array{Any}(undef, np)
    @sync begin
        for i=1:np
            subfreqTxIdx = freqTxIdx[fwdIdx[i]]
            @async epMat[i] = remotecall_fetch(compPrimaryEfieldFreqTx, pids[i],
                                     emMesh, emData, refModel, subfreqTxIdx)
        end
    end

    # assemble the results
    epField = zeros(ComplexF64, 0)
    for i=1:np
        append!(epField, epMat[i])
    end

    nx, ny, nz = emMesh.gridSize
    nE = nx * (ny+1) * (nz+1) + (nx+1) * ny * (nz+1) + (nx+1) * (ny+1) * nz

    epField = reshape(epField, nE, ns, nFreq)

    return epField

end
