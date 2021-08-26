export readEMData, readCSEMData, readMTData
#-------------------------------------------------------------------------------

"""
Read EM data from file.
"""
function readEMData(datafile::String, probType::String)

    if lowercase(probType) == "mt"
        emData, = readMTData(datafile)
    elseif lowercase(probType) == "csem"
        emData, = readCSEMData(datafile)
    elseif lowercase(probType) == "tem"
        emData, = readTEMData(datafile)
    else
        error("problem type: $(probType) is not defined.")
    end

    return emData

end


#-------------------------------------------------------------------------------
"""
Read CSEM data from file.

Input:
    datafile :: String   - the name of the data file.

Output:
    csemData :: CSEMData
    obsData  :: Vector   - data values.
    dataErr  :: Vector   - data standard errors.

"""
function readCSEMData(datafile::String)

    if isfile(datafile)
        fid = open(datafile, "r")
    else
        error("$(datafile) does not exist, please try again.")
    end

    txLoc = [];    rxLoc = [];   dpLen = []
    txID  = [];    rxID  = [];   freqID = []
    dtID  = [];    obsData = []; dataErr = []
    dataType  = []
    freqArray = []
    phaseCon  = []

    while !eof(fid)

        cline = strip(readline(fid))

        # ignore all comments: empty line, or line preceded with #
        isempty(cline)  && continue
        cline[1] == '#' && continue

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
            ns  = parse(Int, tmp[end])
            nd  = 0
            txLoc = zeros(ns, 5)

            while nd < ns
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
                    rxLoc[nd,j] = parse(Float64, cline[j])
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
            dtID  = zeros(Int, nData)
            freqID = zeros(Int, nData)
            obsData = zeros(ComplexF64, nData)
            dataErr = zeros(nData)

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
                dtID[nd]   = parse(Int, cline[4])
                obsData[nd] = parse(Float64, cline[5]) + parse(Float64, cline[6]) * 1im
                dataErr[nd] = parse(Float64, cline[7])
            end

        end

    end

    close(fid)

    if isempty(dpLen)
        dpLen = 0.0
    end
    if isempty(phaseCon)
        phaseCon = "lead"
    end

    csemData = CSEMData(txLoc, rxLoc, dpLen, phaseCon, freqArray, dataType, txID,
                        rxID, freqID, dtID)

    return csemData, obsData, dataErr

end


#-------------------------------------------------------------------------------
"""
Read MT data from file.

Input:
    datafile :: String   - the name of the data file.

Output:
    mtdata   :: MTData
    obsData  :: Vector   - data values.
    dataErr  :: Vector   - data standard errors.

"""
function readMTData(datafile::String)

    if isfile(datafile)
        fid = open(datafile, "r")
    else
        error("$(datafile) does not exist, please try again.")
    end

    rxLoc = [];    rxID    = [];  freqID = []
    dcID  = [];    obsData = [];  dataErr = []
    dataType = []
    dataComp = []
    freqArray = []
    phaseCon  = []
    isComplex = false
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
            if isComplex
                obsData = zeros(ComplexF64, nData)
            else
                obsData = zeros(Float64, nData)
            end
            dataErr = zeros(nData)

            while nd < nData
                cline = strip(readline(fid))
                while isempty(cline) || cline[1] == '#'
                    cline = strip(readline(fid))
                end
                cline = split(cline)
                nd = nd + 1
                freqID[nd] = parse(Int, cline[1])
                rxID[nd]   = parse(Int, cline[2])
                dcID[nd]   = parse(Int, cline[3])
                if isComplex
                    obsData[nd] = parse(Float64, cline[4]) + parse(Float64, cline[5]) * 1im
                    dataErr[nd] = parse(Float64, cline[6])
                else
                    obsData[nd] = parse(Float64, cline[4])
                    dataErr[nd] = parse(Float64, cline[5])
                end

            end

        end

    end

    close(fid)

    if isempty(phaseCon)
        phaseCon = "lead"
    end

    mtdata = MTData(rxLoc, phaseCon, freqArray, dataType, dataComp, rxID, freqID, dcID)

    return mtdata, obsData, dataErr

end
