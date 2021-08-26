
function writeMT1DAniResp(respfile::String, Z::Array, Rho::Array, Phs::Array)

    fid = open(respfile, "w")

    @printf(fid,"%-18s %s\n","# Format:","MT1DAniResp_1.0")
    descrb  = "Response file generated at " * Libc.strftime(time())
    @printf(fid,"%-18s %s\n","# Description:", descrb)

    # frequencies
    nFreq = length(freqs)
    @printf(fid,"%-20s %4d\n","Frequencies (Hz):", nFreq)
    for i = 1:nFreq
        @printf(fid,"%15.5e\n", freqs[i])
    end

    # Impedance
    @printf(fid, "Impedance:\n")
    @printf(fid,"#  %-15s %-15s %-15s %-15s %-15s %-15s %-15s %s\n",
            "RealZXX","ImagZXX", "RealZXY","ImagZXY","RealZYX","ImagZYX","RealZYY","ImagZYY")

    for j in 1:nFreq
        @printf(fid,"%15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e\n",
                real(Z[1, j]), imag(Z[1, j]), real(Z[2, j]), imag(Z[2, j]),
                real(Z[3, j]), imag(Z[3, j]), real(Z[4, j]), imag(Z[4, j]))
    end


    # Apparent resistivity and phase
    @printf(fid, "Apparent resistivity:\n")
    @printf(fid,"#       %-13s %-13s %-13s %-13s %-13s %-13s %-13s %s\n",
            "RhoXX","PhsXX","RhoXY","PhsXY","RhoYX","PhsYX","RhoYY","PhsYY")

    for j in 1:nFreq
        @printf(fid,"%13.4f %13.2f %13.4f %13.2f %13.4f %13.2f %13.4f %13.2f\n",
                Rho[1, j], Phs[1, j], Rho[2, j], Phs[2, j], Rho[3, j], Phs[3, j],
                Rho[4, j], Phs[4, j])
    end

    close(fid)

end
