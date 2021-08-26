export CSEMDataInfo, EMfieldResp, readCSEMResp

mutable struct CSEMDataInfo{T<:Real}

    txLoc     :: Array{T}                   # transmitter location array
    rxLoc     :: Array{T}                   # receiver location array
    phaseCon  :: String                     # phase convention
    freqArray :: Vector{T}                  # frequency array
    dataType  :: Vector{String}             # data type
    txID      :: Vector{Int}                # transmitter index
    rxID      :: Vector{Int}                # receiver index
    freqID    :: Vector{Int}                # frequency index

end


mutable struct EMfieldResp{T<:Complex}

    Ex  :: Array{T}
    Ey  :: Array{T}
    Ez  :: Array{T}
    Bx  :: Array{T}
    By  :: Array{T}
    Bz  :: Array{T}

end

#-------------------------------------------------------------------------------
"""
Read CSEM forward responses from file.

Input:
    respfile :: String   - the name of the response file.

Output:
    dataInfo :: CSENDataInfo
    resp     :: EMfieldResp   - forward response.

"""
function readCSEMResp(respfile::String)

    if isfile(respfile)
        fid = open(respfile, "r")
    else
        error("$(respfile) does not exist, please try again.")
    end

    txLoc = [];    rxLoc = [];   dpLen = []
    txID  = [];    rxID  = [];   freqID = []
    dataType  = []
    freqArray = []
    phaseCon  = []

    nFreq = [];  nTx = [];  nRx = []

    # Ex = [];  Ey = [];  Ez = []
    # Bx = [];  By = [];  Bz = []
    respTmp = []

    while !eof(fid)

        cline = strip(readline(fid))

        # ignore all comments: empty line, or line preceded with #
        while isempty(cline) || cline[1] == '#'
            cline = strip(readline(fid))
        end

        # data format
        if occursin("Format", cline)
            tmp = split(cline)
            format = tmp[2]

        # source type
        elseif occursin("Source Type", cline)
            tmp = split(cline)
            srcType = tmp[end]

        # dipole length, optional
        elseif occursin("Dipole Length", cline)
            tmp = split(cline)
            dpLen = float(tmp[end])

        # phase convention
        elseif occursin("Phase Convention", cline)
            tmp = split(cline)
            phaseCon = string(tmp[end])

        # source location
        elseif occursin("Source Location", cline)
            tmp = split(cline)
            nTx  = parse(Int, tmp[end])
            nd  = 0
            txLoc = zeros(nTx, 5)

            while nd < nTx
                cline = strip(readline(fid))
                while isempty(cline) || cline[1] == '#'
                    cline = strip(readline(fid))
                end
                cline = split(cline)
                nd = nd + 1
                for j = 1:5
                    txLoc[nd,j] = parse(Float64, cline[j])
                end
            end

        # receiver location
        elseif occursin("Receiver Location", cline)
            tmp = split(cline)
            nRx  = parse(Int, tmp[end])
            nd  = 0
            rxLoc = zeros(nRx, 3)
            while nd < nRx
                cline = strip(readline(fid))
                while isempty(cline) || cline[1] == '#'
                    cline = strip(readline(fid))
                end
                cline = split(cline)
                nd = nd + 1
                for j = 1:3
                    rxLoc[nd, j] = parse(Float64, cline[j])
                end
            end

        # frequencies
        elseif occursin("Frequencies", cline)
            tmp = split(cline)
            nFreq  = parse(Int, tmp[end])
            nd  = 0
            freqArray = zeros(nFreq)
            while nd < nFreq
                cline = strip(readline(fid))
                while isempty(cline) || cline[1] == '#'
                    cline = strip(readline(fid))
                end
                #cline = split(cline)
                nd = nd + 1
                freqArray[nd] = parse(Float64, cline)
            end

        # data type
        elseif occursin("DataType", cline)
            tmp = split(cline)
            nDt = parse(Int, tmp[end])
            dataType = Array{String}(undef, nDt)
            nd = 0
            while nd < nDt
                cline = strip(readline(fid))
                nd = nd + 1
                dataType[nd] = cline
            end

        # data block
        elseif occursin("Data Block", cline)
            tmp   = split(cline)
            nData = parse(Int, tmp[end])
            nd    = 0
            txID  = zeros(Int, nData)
            rxID  = zeros(Int, nData)
            freqID = zeros(Int, nData)
            respTmp = zeros(ComplexF64, nData, 6)

            while nd < nData
                cline = strip(readline(fid))
                while isempty(cline) || cline[1] == '#'
                    cline = strip(readline(fid))
                end
                cline = split(cline)
                nd = nd + 1
                freqID[nd] = parse(Int, cline[1])
                txID[nd]   = parse(Int, cline[2])
                rxID[nd]   = parse(Int, cline[3])

                for j=1:6
                    respTmp[nd, j] = parse(Float64, cline[2*j+2]) + parse(Float64, cline[2*j+3]) * 1im
                end

            end

        end

    end

    close(fid)

    if isempty(phaseCon)
        phaseCon = "lead"
    end

    Ex = reshape(respTmp[:, 1], nRx, nTx, nFreq)
    Ey = reshape(respTmp[:, 2], nRx, nTx, nFreq)
    Ez = reshape(respTmp[:, 3], nRx, nTx, nFreq)
    Bx = reshape(respTmp[:, 4], nRx, nTx, nFreq)
    By = reshape(respTmp[:, 5], nRx, nTx, nFreq)
    Bz = reshape(respTmp[:, 6], nRx, nTx, nFreq)

    dataInfo = CSEMDataInfo(txLoc, rxLoc, phaseCon, freqArray, dataType, txID, rxID, freqID)
    resp     = EMfieldResp(Ex, Ey, Ez, Bx, By, Bz)

    return dataInfo, resp

end
