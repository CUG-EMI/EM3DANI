export writeEMFwdResp
#-------------------------------------------------------------------------------
"""
Write CSEM forward responses into file.

Input:
    respfile :: String   - the name of the response file for output.
    emdata   :: CSEMData
    fwdResp  :: Array    - forward response.

"""
function writeEMFwdResp(respfile::String, emData::CSEMData, fwdResp::Array)

    fid = open(respfile, "w")

    @printf(fid,"%-18s %s\n","# Format:","CSEMResp_1.0")
    descrb  = "Response file generated at " * Libc.strftime(time())
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


    nData  = nFreq * nTx * nRx

    @printf(fid,"%-15s %d\n","Data Block:",nData)


    @printf(fid,"# %-8s %-7s %-17s %-29s %-29s %-29s %-29s %-29s %s\n","FreqNo.","TxNo.","RxNo.",
            "Ex(Re,Im)","Ey(Re,Im)","Ez(Re,Im)","Bx(Re,Im)","By(Re,Im)","Bz(Re,Im)")

    for i = 1:nFreq, j=1:nTx, k=1:nRx
        @printf(fid,"%6d %8d %7d %17.6e %13.6e %15.6e %13.6e %15.6e %13.6e %15.6e %13.6e %15.6e %13.6e %15.6e %13.6e\n",
                i, j, k, real(fwdResp[k,j,i]), imag(fwdResp[k,j,i]), real(fwdResp[k,nTx+j,i]), imag(fwdResp[k,nTx+j,i]),
                real(fwdResp[k,2*nTx+j,i]), imag(fwdResp[k,2*nTx+j,i]), real(fwdResp[k,3*nTx+j,i]), imag(fwdResp[k,3*nTx+j,i]),
                real(fwdResp[k,4*nTx+j,i]), imag(fwdResp[k,4*nTx+j,i]), real(fwdResp[k,5*nTx+j,i]), imag(fwdResp[k,5*nTx+j,i]) )
    end

    close(fid)

end



"""
Write MT forward responses into file.

Input:
    respfile :: String   - the name of the response file for output.
    emdata   :: MTData
    fwdResp  :: Array    - forward response.

"""
function writeEMFwdResp(respfile::String, emData::MTData, fwdResp::Array)

    fid = open(respfile, "w")

    @printf(fid,"%-18s %s\n","# Format:","MT3DResp_1.0")
    descrb  = "Response file generated at " * Libc.strftime(time())
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


    nData  = nFreq * nRx

    @printf(fid,"%-15s %d\n","Data Block:",nData)


    if dataType == "Impedance"
        @printf(fid,"# %-8s %-9s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %s\n","FreqNo.","RxNo.",
                "RealZXX","ImagZXX", "RealZXY","ImagZXY","RealZYX","ImagZYX","RealZYY","ImagZYY")

        for j in 1:nFreq, k in 1:nRx
            @printf(fid,"%6d %8d %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e\n", j, k,
                    real(fwdResp[k,1,j]), imag(fwdResp[k,1,j]), real(fwdResp[k,2,j]), imag(fwdResp[k,2,j]),
                    real(fwdResp[k,3,j]), imag(fwdResp[k,3,j]), real(fwdResp[k,4,j]), imag(fwdResp[k,4,j]))
        end

    elseif dataType == "Impedance_Tipper"
        @printf(fid,"# %-8s %-9s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %s\n",
                "FreqNo.","RxNo.","RealZXX","ImagZXX", "RealZXY","ImagZXY","RealZYX","ImagZYX",
                "RealZYY","ImagZYY","RealTZX","ImagTZX","RealTZY","ImagTZY")

        for j in 1:nFreq, k in 1:nRx
            @printf(fid,"%6d %8d %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e\n", j, k,
                    real(fwdResp[k,1,j]), imag(fwdResp[k,1,j]), real(fwdResp[k,2,j]), imag(fwdResp[k,2,j]),
                    real(fwdResp[k,3,j]), imag(fwdResp[k,3,j]), real(fwdResp[k,4,j]), imag(fwdResp[k,4,j]),
                    real(fwdResp[k,5,j]), imag(fwdResp[k,5,j]), real(fwdResp[k,6,j]), imag(fwdResp[k,6,j]))
        end

    elseif dataType == "Rho_Phs"
        @printf(fid,"# %-8s %-12s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %s\n","FreqNo.","RxNo.",
                "RhoXX","PhsXX","RhoXY","PhsXY","RhoYX","PhsYX","RhoYY","PhsYY")

        for j in 1:nFreq, k in 1:nRx
            @printf(fid,"%6d %8d %13.2f %13.2f %13.2f %13.2f %13.2f %13.2f %13.2f %13.2f\n", j, k,
                    fwdResp[k,1,j], fwdResp[k,2,j], fwdResp[k,3,j], fwdResp[k,4,j],
                    fwdResp[k,5,j], fwdResp[k,6,j], fwdResp[k,7,j], fwdResp[k,8,j])
        end

    elseif dataType == "Rho_Phs_Tipper"
        @printf(fid,"# %-8s %-12s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-15s %-15s %-15s %s\n","FreqNo.","RxNo.",
                "RhoXX","PhsXX","RhoXY","PhsXY","RhoYX","PhsYX","RhoYY","PhsYY",
                "RealTZX","ImagTZX","RealTZY","ImagTZY")

        for j in 1:nFreq, k in 1:nRx
            @printf(fid,"%6d %8d %13.2f %13.2f %13.2f %13.2f %13.2f %13.2f %13.2f %13.2f %15.6e %15.6e %15.6e %15.6e\n", j, k,
                    fwdResp[k,1,j], fwdResp[k,2,j], fwdResp[k,3,j], fwdResp[k,4,j],
                    fwdResp[k,5,j], fwdResp[k,6,j], fwdResp[k,7,j], fwdResp[k,8,j],
                    fwdResp[k,9,j], fwdResp[k,10,j], fwdResp[k,11,j], fwdResp[k,12,j])
        end

    end  # dataType==?

    close(fid)

end
