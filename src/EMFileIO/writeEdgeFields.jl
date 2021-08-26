export writeEdgeFields
#-------------------------------------------------------------------------------

"""
Write grid edge vector variables (e.g. E-field) into file.

Input:
    fieldfile :: String   - the name of the file for output.
    eField    :: Array    - grid edge vector variables.
    gridSize  :: Vector   - cell numbers in x,y,z direction.

"""
function writeEdgeFields(fieldfile::String, eField::Array{ComplexF64},
                         gridSize::Vector{Int})

    nE    = size(eField, 1)
    nTx   = size(eField, 2)
    nFreq = size(eField, 3)

    nx, ny, nz = gridSize
    nEx = nx * (ny+1) * (nz+1)
    nEy = (nx+1) * ny * (nz+1)
    nEz = (nx+1) * (ny+1) * nz

    if nE != nEx+nEy+nEz
        error("Grid dimension dismatch in function writeEdgeFields!")
    end

    fid = open(fieldfile, "w")
    @printf(fid,"%s %5d\n","# Frequencies:", nFreq)
    @printf(fid,"%s %5d\n","# Transmitters:", nTx)
    @printf(fid,"%s %5d %5d %5d\n","Grid dimension:", nx,ny,nz)
    @printf(fid,"\n")

    for i=1:nFreq
        @printf(fid,"%s %4d\n","Frequency #:", i)

        for j=1:nTx
            @printf(fid,"%s %4d\n","Transmitter #:", j)

            # x-component
            @printf(fid,"x-component: real\n")
            for k=1:nEx
                @printf(fid,"%15.6e", real(eField[k,j,i]))
                if mod(k,20) == 0; @printf(fid,"\n"); end
            end
            if mod(nEx,20) > 0; @printf(fid,"\n"); end

            @printf(fid,"x-component: imag\n")
            for k=1:nEx
                @printf(fid,"%15.6e", imag(eField[k,j,i]))
                if mod(k,20) == 0; @printf(fid,"\n"); end
            end
            if mod(nEx,20) > 0; @printf(fid,"\n"); end

            # y-component
            @printf(fid,"y-component: real\n")
            for k=1:nEy
                @printf(fid,"%15.6e", real(eField[k+nEx,j,i]))
                if mod(k,20) == 0; @printf(fid,"\n"); end
            end
            if mod(nEy,20) > 0; @printf(fid,"\n"); end

            @printf(fid,"y-component: imag\n")
            for k=1:nEy
                @printf(fid,"%15.6e", imag(eField[k+nEx,j,i]))
                if mod(k,20) == 0; @printf(fid,"\n"); end
            end
            if mod(nEy,20) > 0; @printf(fid,"\n"); end

            # z-component
            @printf(fid,"z-component: real\n")
            for k=1:nEz
                @printf(fid,"%15.6e", real(eField[k+nEx+nEy,j,i]))
                if mod(k,20) == 0; @printf(fid,"\n"); end
            end
            if mod(nEz,20) > 0; @printf(fid,"\n"); end

            @printf(fid,"z-component: imag\n")
            for k=1:nEz
                @printf(fid,"%15.6e", imag(eField[k+nEx+nEy,j,i]))
                if mod(k,20) == 0; @printf(fid,"\n"); end
            end
            if mod(nEz,20) > 0; @printf(fid,"\n"); end

        end   # over nTx
    end   # over nFreq

    close(fid)

end
