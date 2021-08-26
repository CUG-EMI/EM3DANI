export writeCSEMData, writeMTData
#-------------------------------------------------------------------------------
"""
Write predicted CSEM data into file.

Input:
    datafile :: String   - the name of the data file for output.
    emdata   :: CSEMData
    predData :: Vector   - data values.
    dataErr  :: Vector   - data standard errors (optional).

"""
function writeCSEMData(datafile::String, emData::CSEMData, predData::Vector,
                       dataErr::Vector{Float64}=zeros(0))

    fid = open(datafile, "w")

    # data format
    @printf(fid,"%-18s%s\n", "# Format:", "CSEM3DData_1.0")
    descrb  = "Data file generated at " * Libc.strftime(time())
    @printf(fid,"%-18s %s\n","# Description:", descrb)

    # source type
    srcType = "HED"
    @printf(fid,"%-20s %s\n","Source Type:",srcType)

    # dipole length
    if emData.dpLen > 1e-1
        @printf(fid,"%-20s %G\n","Dipole Length:", emData.dpLen)
    end

    # phase convention
    phaseCon = emData.phaseCon
    @printf(fid,"%-20s %s\n", "Phase Convention:", phaseCon)

    # source location
    txLoc = emData.txLoc
    nTx = size(txLoc, 1)
    @printf(fid,"%-25s %4d\n","Source Location (m):", nTx)
    @printf(fid,"%-10s %-12s %-12s %-6s %-16s %s\n","#","X","Y","Z","Azimuth","Dip")
    for i = 1:nTx
        for j = 1:5
            @printf(fid,"%12.2f ",txLoc[i,j])
        end
        @printf(fid,"\n")
    end

    # receiver location
    rxLoc = emData.rxLoc
    nRx = size(rxLoc, 1)
    @printf(fid,"%-25s %4d\n","Receiver Location (m):", nRx)
    @printf(fid,"%-10s %-12s %-12s %s\n","#","X","Y","Z")
    for i = 1:nRx
        @printf(fid,"%12.2f %12.2f %12.2f\n", rxLoc[i,1], rxLoc[i,2], rxLoc[i,3])
    end

    # frequencies
    freqs = emData.freqArray
    nFreq = length(freqs)
    @printf(fid,"%-20s %4d\n","Frequencies (Hz):", nFreq)
    for i = 1:nFreq
        @printf(fid,"%15.5e\n", freqs[i])
    end

    # data type
    dataType = emData.dataType
    nDT = length(dataType)
    @printf(fid,"DataType:  %4d\n", nDT)
    for i = 1:nDT
        @printf(fid,"%s\n",dataType[i])
    end

    # data block
    # DataType includes:
    #            Ex (real and imag), Ey(real and imag), Ez(real and imag)
    #            Bx (real and imag), By(real and imag), Bz(real and imag)
    if isempty(dataErr)
        dataErr = abs.(predData) * 0.02
    elseif length(dataErr) .== 1
        dataErr = abs.(predData) * dataErr[1]
    end

    txID   = emData.txID
    rxID   = emData.rxID
    freqID = emData.freqID
    dtID   = emData.dtID
    nData  = size(predData,1)

    @printf(fid,"%-15s %d\n", "Data Block:", nData)
    # freqID   scrID   recID   DataTpye  RealPart ImagPart DataError
    @printf(fid,"# %-8s %-7s %-7s %-9s %-15s %-15s %s\n","FreqNo.","TxNo.","RxNo.",
            "DTypeNo.","RealValue","ImagValue","Error")
    for i = 1:nData
        @printf(fid,"%6d %8d %7d %8d %15.6e %15.6e %15.6e\n", freqID[i],
                txID[i], rxID[i], dtID[i], real(predData[i]), imag(predData[i]), dataErr[i])
    end

    close(fid)

end


#-------------------------------------------------------------------------------
"""
Write predicted MT data into file.

Input:
    datafile :: String   - the name of the data file for output.
    emdata   :: MTData
    predData :: Vector   - data values.
    dataErr  :: Vector   - data standard errors (optional).

"""
function writeMTData(datafile::String, emData::MTData, predData::Array, dataErr=zeros(0))

    fid = open(datafile, "w")

    # data format
    @printf(fid, "%-18s%s\n", "# Format:", "MT3DData_1.0")
    descrb  = "Data file generated at " * Libc.strftime(time())
    @printf(fid,"%-18s %s\n","# Description:", descrb)

    # phase convention
    phaseCon = emData.phaseCon
    @printf(fid,"%-20s %s\n", "Phase Convention:", phaseCon)

    # receiver location
    rxLoc = emData.rxLoc
    nRx = size(rxLoc, 1)
    @printf(fid,"%-25s %4d\n","Receiver Location (m):", nRx)
    @printf(fid,"%-10s %-12s %-12s %s\n","#","X","Y","Z")
    for i = 1:nRx
        @printf(fid,"%12.2f %12.2f %12.2f\n", rxLoc[i,1], rxLoc[i,2], rxLoc[i,3])
    end

    # frequencies
    freqs = emData.freqArray
    nFreq = length(freqs)
    @printf(fid,"%-20s %4d\n","Frequencies (Hz):", nFreq)
    for i = 1:nFreq
        @printf(fid,"%15.5e\n", freqs[i])
    end

    # data type, data components
    dataType = emData.dataType
    dataComp = emData.dataComp
    nDC = length(dataComp)
    @printf(fid,"DataType:  %s\n", dataType)
    @printf(fid,"DataComp:  %4d\n", nDC)
    for i = 1:nDC
        @printf(fid,"%s\n",dataComp[i])
    end

    # data block
    if isempty(dataErr)
        dataErr = abs.(predData) * 0.05
    elseif length(dataErr) .== 1
        dataErr = abs.(predData) * dataErr[1]
    end

    rxID   = emData.rxID
    freqID = emData.freqID
    dtID   = emData.dtID
    nData  = size(predData, 1)

    @printf(fid,"%-15s %d\n","Data Block:", nData)

    if eltype(predData) <: Complex     # iseltype(predData, Complex)
        @printf(fid,"# %-8s %-7s %-9s %-15s %-15s %s\n","FreqNo.","RxNo.",
                "DCompNo.","RealValue","ImagValue","Error")
        for i = 1:nData
            @printf(fid,"%6d %8d %7d %15.6e %15.6e %15.6e\n", freqID[i], rxID[i],
                    dcID[i], real(predData[i]), imag(predData[i]), dataErr[i])
        end
    else
        @printf(fid,"# %-8s %-7s %-11s %-14s %s\n","FreqNo.","RxNo.",
                "DCompNo.","Value","Error")
        for i = 1:nData
            @printf(fid,"%6d %8d %7d %15.6e %15.6e\n", freqID[i], rxID[i],
                    dcID[i], predData[i], dataErr[i])
        end
    end

    close(fid)

end
