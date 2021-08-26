export initRefModel, get3dBGModel
export compPrimaryElectricField, getPrimaryFieldsAtRxs

#-------------------------------------------------------------------------------
"""
    initRefModel(emMesh; sigEarth, sigAir, sigWater)

Obtain 1D background model for primary filed computation of CSEM.

# Arguments
 - Input
    * `emMesh`   =::CSEMTensorMesh
    * `refModel` =::RefModel, reference model

 - Output
    * `sig1D`   =::Vector, conductivities of background model
    * `depth1D` =::Vector, location of top of each layer

"""
function initRefModel(emMesh::EMTensorMesh; sigAir=1e-8, sigWater=3.3, sigEarth=1.0)

    #
    zLen = emMesh.zLen

    # background model
    sig1D   = zeros(Float64, 0)
    depth1D = zeros(Float64, 0) # location of top of layer

    # locate boundary of background model
    # check if air layer is present
    nAir = 0
    if !isempty(emMesh.airLayer)

        nAir = length(emMesh.airLayer)
        airDepth = sum(emMesh.airLayer)
        topDepth = -emMesh.origin[3]
        append!(sig1D, [sigAir])
        append!(depth1D, [topDepth])

    end

    # check if seawater is present
    nSea = 0
    if !isempty(emMesh.seaLayer)

        nSea = length(emMesh.seaLayer)
        seaDepth = sum(emMesh.seaLayer)
        append!(sig1D, [sigWater])
        if nAir > 0
            topDepth = airDepth - emMesh.origin[3]
        else
            topDepth = -emMesh.origin[3]
        end
        append!(depth1D, [topDepth])
    end

    # uniform sediment
    append!(sig1D, collect(sigEarth))
    topDepth = sum(zLen[1:nAir+nSea]) - emMesh.origin[3]
    append!(depth1D, collect(topDepth))

    return RefModel(sig1D, depth1D)

end


"""
    get3dBGModel(emMesh, refModel)

generates discretized background model for CSEM computation.

"""
function get3dBGModel(emMesh::EMTensorMesh, refModel::RefModel)

    zNode = cumsum([0; emMesh.zLen]) .- emMesh.origin[3]

    sig1D   = copy(refModel.sig1D)
    depth1D = copy(refModel.depth1D)

    nx, ny, nz = emMesh.gridSize
    sig3DBG = ones(nx, ny, nz) * sig1D[end]

    for j=1:length(depth1D)-1
        idx1 = findfirst(isequal(depth1D[j]), zNode)
        idx2 = findfirst(isequal(depth1D[j+1]), zNode)
        sig3DBG[:, :, idx1:idx2-1] .= sig1D[j]
    end

    #
    sig3DBG = vec(sig3DBG)
    if emMesh.condType == "anisotropy"
        sig3DBG = hcat(sig3DBG, sig3DBG, sig3DBG)
    end

    return sig3DBG

end


"""
    compPrimaryElectricField(emMesh, txLoc, freq, sigEarth, sigAir, sigWater)

calculates quasi-analytic primary electric field for layered model.

# Arguments
 - Input
    * `emMesh`   =::EMTensorMesh
    * `txLoc`    =::Array[1x3], transmitter location
    * `freq`     =::Array, frequencies
    * `refModel` =::RefModel, reference model

 - Output
    * `e1d`  =::Array[nE x nFreq], primary electric field

"""
function compPrimaryElectricField(emMesh::EMTensorMesh, txLoc::Array{T},
            freq::Array{T}, dipLen::T, refModel::RefModel) where {T<:Real}

    # grid location
    xNode = cumsum([0; emMesh.xLen]) .- emMesh.origin[1]
    yNode = cumsum([0; emMesh.yLen]) .- emMesh.origin[2]
    zNode = cumsum([0; emMesh.zLen]) .- emMesh.origin[3]

    # get background model
    sig1D   = copy(refModel.sig1D)
    depth1D = copy(refModel.depth1D)

    # get sampling locations of EM fields
    # first get cell center coordiante
    xCen = xNode[1:end-1] + diff(xNode) / 2
    yCen = yNode[1:end-1] + diff(yNode) / 2
    zCen = zNode[1:end-1] + diff(zNode) / 2

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

    nE = nx * (ny+1) * (nz+1) + (nx+1) * ny * (nz+1) + (nx+1) * (ny+1) * nz
    ns = size(txLoc, 1)
    nFreq = length(freq)
    nez   = size(ezLoc, 1)
    e1d   = zeros(ComplexF64, nE, ns, nFreq)

    txtmp  = zeros(1, 5)
    rxtmpz = zeros(1, 3)
    for i = 1:nFreq
        tmpF = collect(freq[i])
        @inbounds for j = 1:ns
            txtmp[1, :] = txLoc[j, :]

            println("   - Freqs. No.: $i / $nFreq,  Tx No.: $j / $ns")
            if dipLen < 10.0
                compFlag = 1
                @printf("%-12s", "     - Ex:")
                t = @elapsed ex1d = comp1DCSEM(txtmp, exLoc, tmpF, sig1D, depth1D, compFlag)
                @printf("%8.2f %s\n", t, "seconds.")
                ex1d = vcat(pEx, vec(ex1d))

                compFlag = 2
                @printf("%-12s", "     - Ey:")
                t = @elapsed ey1d = comp1DCSEM(txtmp, eyLoc, tmpF, sig1D, depth1D, compFlag)
                @printf("%8.2f %s\n", t, "seconds.")
                ey1d = vcat(pEy, vec(ey1d))

                compFlag = 3
                @printf("%-12s", "     - Ez:")
                ez1d = zeros(ComplexF64, 0)
                t = @elapsed begin
                    @inbounds for k = 1:nez
                        rxtmpz[1, :] = ezLoc[k, :]
                        tmp = comp1DCSEM(txtmp, rxtmpz, tmpF, sig1D, depth1D, compFlag)
                        append!(ez1d, collect(tmp))
                    end

                end
                @printf("%8.2f %s\n", t, "seconds.")
                ez1d = vcat(pEz, vec(ez1d))

            else # finite length
                println("Finite length source is used!")
                nIntPts = 15
                compFlag = 1
                @printf("%-12s", "     - Ex:")
                t = @elapsed ex1d = dipLen * comp1DCSEM(txtmp, exLoc, dipLen,
                                       nIntPts, tmpF, sig1D, depth1D, compFlag)
                @printf("%8.2f %s\n", t, "seconds.")
                ex1d = vcat(pEx, vec(ex1d))

                compFlag = 2
                @printf("%-12s", "     - Ey:")
                t = @elapsed ey1d = dipLen * comp1DCSEM(txtmp, eyLoc, dipLen,
                                       nIntPts, tmpF, sig1D, depth1D, compFlag)
                @printf("%8.2f %s\n", t, "seconds.")
                ey1d = vcat(pEy, vec(ey1d))

                compFlag = 3
                @printf("%-12s", "     - Ez:")
                ez1d = zeros(ComplexF64, 0)
                t = @elapsed begin
                    @inbounds for k = 1:nez
                        rxtmpz[1, :] = ezLoc[k, :]
                        tmp = dipLen * comp1DCSEM(txtmp, rxtmpz,
                              dipLen, nIntPts, tmpF, sig1D, depth1D, compFlag)
                        append!(ez1d, collect(tmp))
                    end
                end
                @printf("%8.2f %s\n", t, "seconds.")
                ez1d = vcat(pEz, vec(ez1d))
            end
            e1d[:, j, i] = vcat(ex1d, ey1d, ez1d)
        end
    end

    return e1d

end


"""
    compPrimaryElectricField(emMesh, emData, refModel)

calculates quasi-analytic primary electric field for layered model.

"""
function compPrimaryElectricField(emMesh::EMTensorMesh, emData::CSEMData,
                                  refModel::RefModel, freqIdx=zeros(0))

    #
    txLoc = emData.txLoc
    dipLen = emData.dpLen
    if isempty(freqIdx)
        freqArray = emData.freqArray
    else
        freqArray = emData.freqArray[freqIdx]
    end

    epField = compPrimaryElectricField(emMesh, txLoc, freqArray, dipLen, refModel)

    return epField

end


#
# This routine is particularly for MT, while the above one is for CSEM.
#
function compPrimaryElectricField(emMesh::EMTensorMesh, emData::MTData,
                                  refModel::RefModel, freqIdx=zeros(0))

    #
    sig1D    = refModel.sig1D
    depth1D  = refModel.depth1D
    nx, ny, nz = emMesh.gridSize

    if isempty(freqIdx)
        freqArray = emData.freqArray
    else
        freqArray = emData.freqArray[freqIdx]
    end

    zNode = cumsum([0; emMesh.zLen]) .- emMesh.origin[3]

    sigLayer = ones(nz) * sig1D[end]
    for j=1:length(depth1D)-1
        idx1 = findfirst(isequal(depth1D[j]), zNode)
        idx2 = findfirst(isequal(depth1D[j+1]), zNode)
        sigLayer[idx1:idx2-1] .= sig1D[j]
    end

    nFreq = length(freqArray)
    nEx = nx * (ny+1) * (nz+1)
    nEy = (nx+1) * ny * (nz+1)
    nEz = (nx+1) * (ny+1) * nz
    nE = nEx + nEy + nEz

    epField = zeros(ComplexF64, nE, 2, nFreq)

    for iFreq=1:nFreq
        Ex = zeros(ComplexF64, nx, ny+1, nz+1)
        Ey = zeros(ComplexF64, nx+1, ny, nz+1)

        epTmp = mt1DAnalyticFields(freqArray[iFreq], sigLayer, zNode)

        for iz = 1:nz+1
            Ex[:, :, iz] .= epTmp[iz]
            Ey[:, :, iz] .= epTmp[iz]
        end
        epField[1:nEx, 1, iFreq] = vec(Ex)
        epField[nEx+1:nEx+nEy, 2, iFreq] = vec(Ey)
    end

    return epField
end


"""
    getPrimaryFieldsAtRxs(emData, refModel)

calculates primary fields at receivers.

#Arguments
 - Input
    * `emData`   =::CSEMData
    * `refModel` =::RefModel

 - Output
    * `predData` =::Array, primary fields at receivers.

"""
function getPrimaryFieldsAtRxs(emData::CSEMData, refModel::RefModel,
                               freqIdx=zeros(0))

    #
    txLoc = copy(emData.txLoc)
    rxLoc = copy(emData.rxLoc)
    if isempty(freqIdx)
        freq = copy(emData.freqArray)
        freqIdx = collect(Int, 1:length(freq))
    else
        freq = emData.freqArray[freqIdx]
    end

    nFreq = length(freq)
    (nFreq == 1) && ( freq = collect(freq) )

    # get background model
    sig1D   = copy(refModel.sig1D)
    depth1D = copy(refModel.depth1D)

    #
    Exp = [];  Eyp = [];  Ezp = []
    Bxp = [];  Byp = [];  Bzp = []
    dipLen = emData.dpLen
    if dipLen < 10.0
        Exp = comp1DCSEM(txLoc, rxLoc, freq, sig1D, depth1D, 1)
        Eyp = comp1DCSEM(txLoc, rxLoc, freq, sig1D, depth1D, 2)
        Ezp = comp1DCSEM(txLoc, rxLoc, freq, sig1D, depth1D, 3)
        Bxp = comp1DCSEM(txLoc, rxLoc, freq, sig1D, depth1D, 4)
        Byp = comp1DCSEM(txLoc, rxLoc, freq, sig1D, depth1D, 5)
        Bzp = comp1DCSEM(txLoc, rxLoc, freq, sig1D, depth1D, 6)

    else
        printstyled("Finite length source is used!\n", color=:cyan)
        nIntPts = 15
        Exp = dipLen * comp1DCSEM(txLoc, rxLoc, dipLen, nIntPts, freq, sig1D, depth1D, 1)
        Eyp = dipLen * comp1DCSEM(txLoc, rxLoc, dipLen, nIntPts, freq, sig1D, depth1D, 2)
        Ezp = dipLen * comp1DCSEM(txLoc, rxLoc, dipLen, nIntPts, freq, sig1D, depth1D, 3)
        Bxp = dipLen * comp1DCSEM(txLoc, rxLoc, dipLen, nIntPts, freq, sig1D, depth1D, 4)
        Byp = dipLen * comp1DCSEM(txLoc, rxLoc, dipLen, nIntPts, freq, sig1D, depth1D, 5)
        Bzp = dipLen * comp1DCSEM(txLoc, rxLoc, dipLen, nIntPts, freq, sig1D, depth1D, 6)
    end

    return Exp, Eyp, Ezp, Bxp, Byp, Bzp

end # getPrimaryFieldsAtRxs


"""
    getPrimaryFieldsAtRxs(emData, refModel, freqIdx, txIdx)

calculates primary fields at receivers.

#Arguments
 - Input
    * `emData`   =::CSEMData
    * `refModel` =::RefModel
    * `freqIdx`  =::Array
    * `txIdx`    =::Array

 - Output
    * `predData` =::Array, primary fields at receivers.

"""
function getPrimaryFieldsAtRxs(emData::CSEMData, refModel::RefModel,
                                 freqIdx::Array{Int}, txIdx::Array{Int})

    #
    txLoc = emData.txLoc[txIdx, :]
    rxLoc = copy(emData.rxLoc)
    freq  = emData.freqArray[freqIdx]

    nFreq = length(freq)
    (nFreq == 1) && ( freq = collect(freq) )

    # get background model
    sig1D   = copy(refModel.sig1D)
    depth1D = copy(refModel.depth1D)

    #
    Exp = [];  Eyp = [];  Ezp = []
    Bxp = [];  Byp = [];  Bzp = []
    dipLen = emData.dpLen
    if dipLen < 10.0
        Exp = comp1DCSEM(txLoc, rxLoc, freq, sig1D, depth1D, 1)
        Eyp = comp1DCSEM(txLoc, rxLoc, freq, sig1D, depth1D, 2)
        Ezp = comp1DCSEM(txLoc, rxLoc, freq, sig1D, depth1D, 3)
        Bxp = comp1DCSEM(txLoc, rxLoc, freq, sig1D, depth1D, 4)
        Byp = comp1DCSEM(txLoc, rxLoc, freq, sig1D, depth1D, 5)
        Bzp = comp1DCSEM(txLoc, rxLoc, freq, sig1D, depth1D, 6)

    else
        printstyled("Finite length source is used!\n", color=:cyan)
        nIntPts = 15
        Exp = dipLen * comp1DCSEM(txLoc, rxLoc, dipLen, nIntPts, freq, sig1D, depth1D, 1)
        Eyp = dipLen * comp1DCSEM(txLoc, rxLoc, dipLen, nIntPts, freq, sig1D, depth1D, 2)
        Ezp = dipLen * comp1DCSEM(txLoc, rxLoc, dipLen, nIntPts, freq, sig1D, depth1D, 3)
        Bxp = dipLen * comp1DCSEM(txLoc, rxLoc, dipLen, nIntPts, freq, sig1D, depth1D, 4)
        Byp = dipLen * comp1DCSEM(txLoc, rxLoc, dipLen, nIntPts, freq, sig1D, depth1D, 5)
        Bzp = dipLen * comp1DCSEM(txLoc, rxLoc, dipLen, nIntPts, freq, sig1D, depth1D, 6)
    end

    return Exp, Eyp, Ezp, Bxp, Byp, Bzp

end # getPrimaryFieldsAtRxs


include("parcompPrimaryField.jl")
