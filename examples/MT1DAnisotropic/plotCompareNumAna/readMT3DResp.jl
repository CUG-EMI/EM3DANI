export MTDataInfo, readMT3DResp

mutable struct MTDataInfo{T<:Real}

    rxLoc     :: Array{T}                   # receiver location array
    phaseCon  :: String                     # phase convention
    freqArray :: Vector{T}                  # frequency array
    dataType  :: String                     # data type
    dataComp  :: Vector{String}             # data component
    rxID      :: Vector{Int}                # receiver index
    freqID    :: Vector{Int}                # frequency index

end

#-------------------------------------------------------------------------------
"""
Read MT forward responses from file.

Input:
    respfile :: String   - the name of the response file.

Output:
    dataInfo :: MTDataInfo
    resp     :: Array   - forward response.

"""
function readMT3DResp(respfile::String)

    if isfile(respfile)
        fid = open(respfile, "r")
    else
        error("$(respfile) does not exist, please try again.")
    end

    rxLoc = [];    rxID    = [];  freqID = []
    dataType = []
    dataComp = []
    freqArray = []
    phaseCon  = []
    isComplex = false

    resp = []

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

        # phase convention
        elseif occursin("Phase Convention", cline)
            tmp = split(cline)
            phaseCon = string(tmp[end])

        # receiver location
        elseif occursin("Receiver Location", cline)
            tmp = split(cline)
            nr  = parse(Int, tmp[end])
            nd  = 0
            rxLoc = zeros(nr, 3)
            while nd < nr
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
            nf  = parse(Int, tmp[end])
            nd  = 0
            freqArray = zeros(nf)
            while nd < nf
                cline = strip(readline(fid))
                while isempty(cline) || cline[1] == '#'
                    cline = strip(readline(fid))
                end
                #cline = split(cline)
                nd = nd + 1
                freqArray[nd] = parse(Float64, cline)
            end

        # data type, allowed values are:
        #    Impedance, Impedance_Tipper, Rho_Rhs, Rho_Rhs_Tipper
        elseif occursin("DataType", cline)
            tmp = split(cline)
            dataType = string(tmp[end])
            if !occursin(dataType, "Impedance") && !occursin(dataType, "Rho_Phs")
                error("$(dataType) is not supported.")
            end

            occursin("Impedance", dataType)  &&  ( isComplex = true )

        # data components, allowed values for each data type:
        # Impedance:         ZXX, ZXY, ZYX, ZYY
        # Impedance_Tipper:  ZXX, ZXY, ZYX, ZYY, TZX, TZY
        # Rho_Phs:           RhoXX, PhsXX, RhoXY, PhsXY, RhoYX, PhsYX, RhoYY, PhsYY
        # Rho_Phs_Tipper:    RhoXX, PhsXX, RhoXY, PhsXY, RhoYX, PhsYX, RhoYY, PhsYY,
        #                    RealTZX, ImagTZX, RealTZY, ImagTZY
        elseif occursin("DataComp", cline)
            tmp = split(cline)
            nDc = parse(Int, tmp[end])
            dataComp = Array{String}(undef, nDc)
            nd = 0
            while nd < nDc
                cline = strip(readline(fid))
                nd = nd + 1
                dataComp[nd] = cline
            end

        # data block
        elseif occursin("Data Block", cline)
            tmp   = split(cline)
            nData = parse(Int, tmp[end])
            nd    = 0
            rxID  = zeros(Int, nData)
            dcID  = zeros(Int, nData)
            freqID = zeros(Int, nData)
            resp = zeros(Float64, nData, 8)

            while nd < nData
                cline = strip(readline(fid))
                while isempty(cline) || cline[1] == '#'
                    cline = strip(readline(fid))
                end
                cline = split(cline)
                nd = nd + 1
                freqID[nd] = parse(Int, cline[1])
                rxID[nd]   = parse(Int, cline[2])

                for j=1:8
                    resp[nd, j] = parse(Float64, cline[j+2])
                end

            end

        end

    end

    close(fid)

    if isempty(phaseCon)
        phaseCon = "lead"
    end

    dataInfo = MTDataInfo(rxLoc, phaseCon, freqArray, dataType, dataComp, rxID, freqID)

    return dataInfo, resp

end
