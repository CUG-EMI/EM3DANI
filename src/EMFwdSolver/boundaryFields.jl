# export getBoundaryMT, getBoundaryMTIso, getBoundaryMTAni, getBoundIndex

"""
Computes boundary E fields by conducting a series of MT 1D forward modeling.

Input:
    freqs    :: Vector - frequency array.
    emMesh   :: EMTensorMesh
    polFlag  :: Int    - indicates polarization mode (1 for Ex-Hy, 2 for Ey-Hx).

Output:
    bc :: Array  - boundary electric fields.
"""
function getBoundaryMT(freqs::Vector{T}, emMesh::EMTensorMesh, polFlag::Int) where {T<:Real}

    if polFlag != 1 && polFlag != 2
        error("Undefined polarization flag: polFlag=$polFlag!")
    end

    if emMesh.condType == "isotropy"
        bc = getBoundaryMTIso(freqs, emMesh, polFlag)
    elseif emMesh.condType == "anisotropy"
        bc = getBoundaryMTAni(freqs, emMesh, polFlag)
    end

    return bc
end


# Isotropic case
function getBoundaryMTIso(freqs::Vector{T}, emMesh::EMTensorMesh, polFlag::Int) where {T<:Real}

    nFreq = length(freqs)

    xLen = emMesh.xLen
    yLen = emMesh.yLen
    zLen = emMesh.zLen
    zNode = [0; cumsum(zLen)]

    nx, ny, nz = emMesh.gridSize

    sigma3D = reshape(emMesh.sigma, nx, ny, nz)

    # number of boundary fields
    nExT = nx * (ny+1)   # nExB = nExT
    nExW = nx * (nz-1)   # nExE = nExW
    nEyT = ny * (nx+1)   # nEyB = nEyT
    nEyS = ny * (nz-1)   # nEyN = nEyS
    nEzS = nz * (ny+1)   # nEzN = nEzS
    nEzW = nz * (nx-1)   # nEzE = nEzW
    nb = 2 * (nExT + nExW + nEyT + nEyS + nEzS + nEzW)

    bc = zeros(ComplexF64, nb, nFreq)


    # for Ex-Hy polarization
    if polFlag == 1

        # top, bottom, west and east boundary E-fields
        # exT = ones(ComplexF64, nx, ny+1, nFreq)
        exB = zeros(ComplexF64, nx, ny+1, nFreq)
        exW = zeros(ComplexF64, nx, nz-1, nFreq)
        exE = copy(exW)

        # 1st compute west and east boundary fields
        for ix=1:nx
            sig1DW = vec(sigma3D[ix, 1, :])
            sig1DE = vec(sigma3D[ix, end, :])

            for j=1:nFreq
                etmp = mt1DAnalyticFields(freqs[j], sig1DW, zNode)
                etmp = etmp / etmp[1]
                exW[ix, :, j] = etmp[2:end-1]
                exB[ix, 1, j] = etmp[end]

                etmp = mt1DAnalyticFields(freqs[j], sig1DE, zNode)
                etmp = etmp / etmp[1]
                exE[ix, :, j]   = etmp[2:end-1]
                exB[ix, end, j] = etmp[end]
            end # nFreq
        end # nx

        # 2nd compute bottom boundary fields
        for iy=2:ny
            yL = yLen[iy-1]
            yR = yLen[iy]

            for ix=1:nx
                sig1DL = vec(sigma3D[ix, iy-1, :])
                sig1DR = vec(sigma3D[ix, iy, :])
                sig1D  = (sig1DL * yL + sig1DR * yR) / (yL + yR)

                for j=1:nFreq
                    etmp = mt1DAnalyticFields(freqs[j], sig1D, zNode)
                    etmp = etmp / etmp[1]
                    exB[ix, iy, j] = etmp[end]
                end # nFreq
            end # nx
        end # ny


        # assign them to bc in order: top, west, east, bottom (be consistent with
        # the order in getBoundIndex)

        # assume electric fields on the top boundary are known (and are constant).
        bc[1:nExT, :] .= 1.0

        for j=1:nFreq
            ks = nExT+1
            ke = nExT+nExW
            bc[ks:ke, j] = vec(exW[:,:,j])

            ks = ke+1
            ke = ke+nExW
            bc[ks:ke, j] = vec(exE[:,:,j])

            ks = ke+1
            ke = ke+nExT
            bc[ks:ke, j] = vec(exB[:,:,j])
        end


    # for Ey-Hx polarization
    elseif polFlag == 2

        # top, bottom, south and north boundary E-fields
        # eyT = ones(ComplexF64, nx+1, ny, nFreq)
        eyB = zeros(ComplexF64, nx+1, ny, nFreq)
        eyS = zeros(ComplexF64, ny, nz-1, nFreq)
        eyN = copy(eyS)

        # 1st compute south and north boundary fields
        for iy=1:ny
            sig1DS = vec(sigma3D[1, iy, :])
            sig1DN = vec(sigma3D[end, iy, :])

            for j=1:nFreq
                etmp = mt1DAnalyticFields(freqs[j], sig1DS, zNode)
                etmp = etmp / etmp[1]
                eyS[iy, :, j] = etmp[2:end-1]
                eyB[1, iy, j] = etmp[end]

                etmp = mt1DAnalyticFields(freqs[j], sig1DN, zNode)
                etmp = etmp / etmp[1]
                eyN[iy, :, j]   = etmp[2:end-1]
                eyB[end, iy, j] = etmp[end]
            end # nFreq
        end # ny

        # 2nd compute bottom boundary fields
        for ix=2:nx
            xL = xLen[ix-1]
            xR = xLen[ix]

            for iy=1:ny
                sig1DL = vec(sigma3D[ix-1, iy, :])
                sig1DR = vec(sigma3D[ix, iy, :])
                sig1D  = (sig1DL * xL + sig1DR * xR) / (xL + xR)

                for j=1:nFreq
                    etmp = mt1DAnalyticFields(freqs[j], sig1D, zNode)
                    etmp = etmp / etmp[1]
                    eyB[ix, iy, j] = etmp[end]
                end # nFreq
            end # ny
        end # nx


        # assign them to bc in order: top, south, north, bottom (be consistent with
        # the order in getBoundIndex)

        # assume electric fields on the top boundary are known (and are constant).
        nbEx = 2 * (nExT + nExW)
        bc[nbEx+1:nbEx+nEyT, :] .= 1.0

        for j=1:nFreq
            ks = nbEx+nEyT+1
            ke = nbEx+nEyT+nEyS
            bc[ks:ke, j] = vec(eyS[:,:,j])

            ks = ke+1
            ke = ke+nEyS
            bc[ks:ke, j] = vec(eyN[:,:,j])

            ks = ke+1
            ke = ke+nEyT
            bc[ks:ke, j] = vec(eyB[:,:,j])
        end

    end # polFlag

    return bc

end



# Anisotropic case
function getBoundaryMTAni(freqs::Vector{T}, emMesh::EMTensorMesh, polFlag::Int) where {T<:Real}

    nFreq = length(freqs)

    xLen = emMesh.xLen
    yLen = emMesh.yLen
    zLen = emMesh.zLen
    zNode = [0; cumsum(zLen)]

    nx, ny, nz = emMesh.gridSize

    sigxx3D = reshape(emMesh.sigma[:, 1], nx, ny, nz)
    sigyy3D = reshape(emMesh.sigma[:, 2], nx, ny, nz)
    sigzz3D = reshape(emMesh.sigma[:, 3], nx, ny, nz)

    if isempty(emMesh.offsigma)
        sigxy3D = zeros(Float64, nx, ny, nz)
        sigxz3D = zeros(Float64, nx, ny, nz)
        sigyz3D = zeros(Float64, nx, ny, nz)
    else
        sigxy3D = reshape(emMesh.offsigma[:, 1], nx, ny, nz)
        sigxz3D = reshape(emMesh.offsigma[:, 2], nx, ny, nz)
        sigyz3D = reshape(emMesh.offsigma[:, 3], nx, ny, nz)
    end


    # number of boundary fields
    nExT = nx * (ny+1)   # nExB = nExT
    nExW = nx * (nz-1)   # nExE = nExW
    nEyT = ny * (nx+1)   # nEyB = nEyT
    nEyS = ny * (nz-1)   # nEyN = nEyS
    nEzS = nz * (ny+1)   # nEzN = nEzS
    nEzW = nz * (nx-1)   # nEzE = nEzW
    nb = 2 * (nExT + nExW + nEyT + nEyS + nEzS + nEzW)

    bc = zeros(ComplexF64, nb, nFreq)


    if polFlag == 1
        eTop = [1.0+0im; 0]     # for Ex-Hy polarization: Ex0=1, Ey0=0

    elseif polFlag == 2
        eTop = [0; 1.0+0im]     # for Ey-Hx polarization: Ex0=0, Ey0=1

    end


    # For west & east boundary, we need to solve Ex and Ez
    exW = zeros(ComplexF64, nx, nz-1, nFreq)
    exE = copy(exW)
    ezW = zeros(ComplexF64, nx-1, nz, nFreq)
    ezE = copy(ezW)

    # For south & north boundary, we need to solve Ey and Ez
    eyS = zeros(ComplexF64, ny, nz-1, nFreq)
    eyN = copy(eyS)
    ezS = zeros(ComplexF64, ny+1, nz, nFreq)
    ezN = copy(ezS)

    # For top & bottom boundary, we need to solve Ex and Ey
    # exT = ones(ComplexF64, nx, ny+1, nFreq)
    # eyT = ones(ComplexF64, nx+1, ny, nFreq)
    exB = zeros(ComplexF64, nx, ny+1, nFreq)
    eyB = zeros(ComplexF64, nx+1, ny, nFreq)


    # will be used to get an averaging Ez
    ezWL = zeros(ComplexF64, nz, nFreq)
    ezWR = copy(ezWL)
    ezEL = copy(ezWL)
    ezER = copy(ezWL)

    ezSL = zeros(ComplexF64, nz, nFreq)
    ezSR = copy(ezSL)
    ezNL = copy(ezSL)
    ezNR = copy(ezSL)

    # Start with W & E boundary, solve Ex and Ez
    for ix=1:nx
        sigxx1D = vec(sigxx3D[ix, 1, :])
        sigyy1D = vec(sigyy3D[ix, 1, :])
        sigzz1D = vec(sigzz3D[ix, 1, :])
        sigxy1D = vec(sigxy3D[ix, 1, :])
        sigxz1D = vec(sigxz3D[ix, 1, :])
        sigyz1D = vec(sigyz3D[ix, 1, :])
        sig1DW = hcat(sigxx1D, sigyy1D, sigzz1D, sigxy1D, sigxz1D, sigyz1D)

        sigxx1D = vec(sigxx3D[ix, end, :])
        sigyy1D = vec(sigyy3D[ix, end, :])
        sigzz1D = vec(sigzz3D[ix, end, :])
        sigxy1D = vec(sigxy3D[ix, end, :])
        sigxz1D = vec(sigxz3D[ix, end, :])
        sigyz1D = vec(sigyz3D[ix, end, :])
        sig1DE = hcat(sigxx1D, sigyy1D, sigzz1D, sigxy1D, sigxz1D, sigyz1D)


        for j=1:nFreq
            extmp, eytmp, eztmpW = mt1DAnalyticFieldsAni(freqs[j], sig1DW, zNode, eTop, timeFac="pos")
            exW[ix, :, j] = extmp[2:end-1]
            exB[ix, 1, j] = extmp[end]

            extmp, eytmp, eztmpE = mt1DAnalyticFieldsAni(freqs[j], sig1DE, zNode, eTop, timeFac="pos")
            exE[ix, :, j]   = extmp[2:end-1]
            exB[ix, end, j] = extmp[end]

            if ix==1
                ezWL[:, j] = eztmpW
                ezEL[:, j] = eztmpE
                continue
            end

            ezWR[:, j] = eztmpW
            ezER[:, j] = eztmpE
        end # nFreq

        if ix==1; continue; end

        # averaging Ez
        wL = xLen[ix] / (xLen[ix-1]+xLen[ix])
        wR = 1-wL

        ezW[ix-1, :, :] = ezWL*wL + ezWR*wR
        ezE[ix-1, :, :] = ezEL*wL + ezER*wR

        ezWL = copy(ezWR)
        ezEL = copy(ezER)

    end   # nx


    # Then deal with S & N boundary, solve Ey and Ez
    for iy=1:ny
        sigxx1D = vec(sigxx3D[1, iy, :])
        sigyy1D = vec(sigyy3D[1, iy, :])
        sigzz1D = vec(sigzz3D[1, iy, :])
        sigxy1D = vec(sigxy3D[1, iy, :])
        sigxz1D = vec(sigxz3D[1, iy, :])
        sigyz1D = vec(sigyz3D[1, iy, :])
        sig1DS = hcat(sigxx1D, sigyy1D, sigzz1D, sigxy1D, sigxz1D, sigyz1D)

        sigxx1D = vec(sigxx3D[end, iy, :])
        sigyy1D = vec(sigyy3D[end, iy, :])
        sigzz1D = vec(sigzz3D[end, iy, :])
        sigxy1D = vec(sigxy3D[end, iy, :])
        sigxz1D = vec(sigxz3D[end, iy, :])
        sigyz1D = vec(sigyz3D[end, iy, :])
        sig1DN = hcat(sigxx1D, sigyy1D, sigzz1D, sigxy1D, sigxz1D, sigyz1D)


        for j=1:nFreq
            extmp, eytmp, eztmpS = mt1DAnalyticFieldsAni(freqs[j], sig1DS, zNode, eTop, timeFac="pos")
            eyS[iy, :, j] = eytmp[2:end-1]
            eyB[1, iy, j] = eytmp[end]

            extmp, eytmp, eztmpN = mt1DAnalyticFieldsAni(freqs[j], sig1DN, zNode, eTop, timeFac="pos")
            eyN[iy, :, j]   = eytmp[2:end-1]
            eyB[end, iy, j] = eytmp[end]

            if iy==1
                ezSL[:, j] = eztmpS
                ezNL[:, j] = eztmpN
                continue
            end

            ezSR[:, j] = eztmpS
            ezNR[:, j] = eztmpN
        end # nFreq

        if iy==1; continue; end

        # averaging Ez
        wL = yLen[iy] / (yLen[iy-1]+yLen[iy])
        wR = 1-wL

        ezS[iy, :, :] = ezSL*wL + ezSR*wR
        ezN[iy, :, :] = ezNL*wL + ezNR*wR

        ezSL = copy(ezSR)
        ezNL = copy(ezNR)

    end   # ny

    ezS[1, :, :]   = ezS[2, :, :]
    ezS[end, :, :] = ezS[end-1, :, :]
    ezN[1, :, :]   = ezN[2, :, :]
    ezN[end, :, :] = ezN[end-1, :, :]


    # Finish with T & B boundary, solve Ex and Ey
    for iy=2:ny
        wL = yLen[iy-1] / (yLen[iy-1]+yLen[iy])
        wR = 1-wL

        for ix=1:nx
            sigxx1D = wL*vec(sigxx3D[ix, iy-1, :]) + wR*vec(sigxx3D[ix, iy, :])
            sigyy1D = wL*vec(sigyy3D[ix, iy-1, :]) + wR*vec(sigyy3D[ix, iy, :])
            sigzz1D = wL*vec(sigzz3D[ix, iy-1, :]) + wR*vec(sigzz3D[ix, iy, :])
            sigxy1D = wL*vec(sigxy3D[ix, iy-1, :]) + wR*vec(sigxy3D[ix, iy, :])
            sigxz1D = wL*vec(sigxz3D[ix, iy-1, :]) + wR*vec(sigxz3D[ix, iy, :])
            sigyz1D = wL*vec(sigyz3D[ix, iy-1, :]) + wR*vec(sigyz3D[ix, iy, :])
            sig1D = hcat(sigxx1D, sigyy1D, sigzz1D, sigxy1D, sigxz1D, sigyz1D)

            for j=1:nFreq
                extmp, = mt1DAnalyticFieldsAni(freqs[j], sig1D, zNode, eTop, timeFac="pos")
                exB[ix, iy, j] = extmp[end]
            end # nFreq
        end # nx
    end # ny

    for ix=2:nx
        wL = xLen[ix-1] / (xLen[ix-1]+xLen[ix])
        wR = 1-wL

        for iy=1:ny
            sigxx1D = wL*vec(sigxx3D[ix-1, iy, :]) + wR*vec(sigxx3D[ix, iy, :])
            sigyy1D = wL*vec(sigyy3D[ix-1, iy, :]) + wR*vec(sigyy3D[ix, iy, :])
            sigzz1D = wL*vec(sigzz3D[ix-1, iy, :]) + wR*vec(sigzz3D[ix, iy, :])
            sigxy1D = wL*vec(sigxy3D[ix-1, iy, :]) + wR*vec(sigxy3D[ix, iy, :])
            sigxz1D = wL*vec(sigxz3D[ix-1, iy, :]) + wR*vec(sigxz3D[ix, iy, :])
            sigyz1D = wL*vec(sigyz3D[ix-1, iy, :]) + wR*vec(sigyz3D[ix, iy, :])
            sig1D = hcat(sigxx1D, sigyy1D, sigzz1D, sigxy1D, sigxz1D, sigyz1D)

            for j=1:nFreq
                extmp, eytmp = mt1DAnalyticFieldsAni(freqs[j], sig1D, zNode, eTop, timeFac="pos")
                eyB[ix, iy, j] = eytmp[end]
            end # nFreq
        end # nx
    end # ny


    # The last step: assign computed boundary fields to 'bc' in order (be consistent with
    # the order in getBoundIndex)

    # Ex: top, west, east, bottom
    bc[1:nExT, :] .= eTop[1]
    for j=1:nFreq
        ks = nExT+1
        ke = nExT+nExW
        bc[ks:ke, j] = vec(exW[:,:,j])

        ks = ke+1
        ke = ke+nExW
        bc[ks:ke, j] = vec(exE[:,:,j])

        ks = ke+1
        ke = ke+nExT
        bc[ks:ke, j] = vec(exB[:,:,j])
    end

    # Ey: top, south, north, bottom
    nbEx = 2 * (nExT + nExW)
    bc[nbEx+1:nbEx+nEyT, :] .= eTop[2]
    for j=1:nFreq
        ks = nbEx+nEyT+1
        ke = nbEx+nEyT+nEyS
        bc[ks:ke, j] = vec(eyS[:,:,j])

        ks = ke+1
        ke = ke+nEyS
        bc[ks:ke, j] = vec(eyN[:,:,j])

        ks = ke+1
        ke = ke+nEyT
        bc[ks:ke, j] = vec(eyB[:,:,j])
    end

    # Ez: south, north, west, east
    nbEh = 2 * (nExT + nExW + nEyT +nEyS)
    for j=1:nFreq
        ks = nbEh+1
        ke = nbEh+nEzS
        bc[ks:ke, j] = vec(ezS[:,:,j])

        ks = ke+1
        ke = ke+nEzS
        bc[ks:ke, j] = vec(ezN[:,:,j])

        ks = ke+1
        ke = ke+nEzW
        bc[ks:ke, j] = vec(ezW[:,:,j])

        ks = ke+1
        ke = ke+nEzW
        bc[ks:ke, j] = vec(ezE[:,:,j])
    end

    return bc

end



"""
`getBoundIndex` splits the indices of edge variables into a boundary part and an
inner part.

Input:
    gridSize :: Vector   - number of girds in three directions.

Output:
    ii :: Vector  - inner index.
    io :: Vector  - boundary index.

"""
function getBoundIndex(gridSize::Vector{Int})

    nx, ny, nz = gridSize
    nEx = nx * (ny+1) * (nz+1)
    nEy = (nx+1) * ny * (nz+1)
    nEz = (nx+1) * (ny+1) * nz

    # for Ex
    idx3D = reshape(collect(1:nEx), nx, ny+1, nz+1)

    iiX = vec(idx3D[:, 2:end-1, 2:end-1])    # inner
    iT = vec(idx3D[:, :, 1])                 # top boundary
    iW = vec(idx3D[:, 1, 2:end-1])           # west boundary
    iE = vec(idx3D[:, end, 2:end-1])         # east boundary
    iB = vec(idx3D[:, :, end])               # bottom boundary
    ioX = vcat(iT, iW, iE, iB)


    # for Ey
    idx3D = reshape(collect(1:nEy), nx+1, ny, nz+1)

    iiY = vec(idx3D[2:end-1, :, 2:end-1])    # inner
    iT = vec(idx3D[:, :, 1])                 # top boundary
    iS = vec(idx3D[1, :, 2:end-1])           # south boundary
    iN = vec(idx3D[end, :, 2:end-1])         # north boundary
    iB = vec(idx3D[:, :, end])               # bottom boundary
    ioY = vcat(iT, iS, iN, iB)
    iiY = iiY .+ nEx
    ioY = ioY .+ nEx


    # for Ez
    idx3D = reshape(collect(1:nEz), nx+1, ny+1, nz)

    iiZ = vec(idx3D[2:end-1, 2:end-1, :])    # inner
    iS = vec(idx3D[1, :, :])                 # south boundary
    iN = vec(idx3D[end, :, :])               # north boundary
    iW = vec(idx3D[2:end-1, 1, :])           # west boundary
    iE = vec(idx3D[2:end-1, end, :])         # east boundary
    ioZ = vcat(iS, iN, iW, iE)
    iiZ = iiZ .+ (nEx + nEy)
    ioZ = ioZ .+ (nEx + nEy)

    ii = vcat(iiX, iiY, iiZ)
    io = vcat(ioX, ioY, ioZ)

    return ii, io

end  # getBoundIndex
