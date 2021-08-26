# export getTotalFieldsAtRxsCSEM, getTotalFieldsAtRxsMT, getTotalFieldsAtRxsMT2
# export compEMRxTerm, getTransFunction, matchData
# export interpFieldsToRx, interpFieldsToRxAni
#-------------------------------------------------------------------------------


# Construct projection matrix for EM receivers.
function compEMRxTerm(rxLoc::Array{T,2}, emMesh::EMTensorMesh) where {T<:Real}


    # receiver projection matrix

    exMat = trilinearInterpEH(1, rxLoc, emMesh)
    eyMat = trilinearInterpEH(2, rxLoc, emMesh)

    #ezMat = trilinearInterpEH(3, rxLoc, emMesh, false)
    #ezMat = trilinearInterpEH(3, rxLoc, emMesh)
    # ezMat = spzeros(0)
    if emMesh.condType == "isotropy"
        ezMat = shantInterpJz(rxLoc, emMesh)
    elseif emMesh.condType == "anisotropy"
        ezMat = trilinearInterpJzAni(rxLoc, emMesh)
    end


    bxMat = trilinearInterpEH(4, rxLoc, emMesh)
    byMat = trilinearInterpEH(5, rxLoc, emMesh)
    bzMat = trilinearInterpEH(6, rxLoc, emMesh)

    MU0 = 4*pi*1e-7
    bxeMat = shantInterpBx(rxLoc, emMesh); bxeMat=MU0*bxeMat
    byeMat = shantInterpBy(rxLoc, emMesh); byeMat=MU0*byeMat

    rxMat = RxProjMat(exMat, eyMat, ezMat, bxMat, byMat, bzMat, bxeMat, byeMat)

    return rxMat

end


"""
Interpolate fields to receiver positions using the receiver projection matrix.

Input:
    emMesh    :: EMTensorMesh
    emData    :: CSEMData
    refModel  :: RefModel
    eField    :: Array   - secondary electric fields at grid sampling points.
    bField    :: Array   - secondary magnetic fields at grid sampling points.
    epField   :: Array   - primary electric fields at grid sampling points.

Output:
    total fields (all six components) at receivers for all freqs, Tx.

"""
function getTotalFieldsAtRxsCSEM(emMesh::EMTensorMesh, emData::CSEMData,
                    refModel::RefModel, eField::Array{T,3}, bField::Array{T,3},
                    args...) where {T<:Complex}


    # deal with optional arguments
    epField = zeros(0)
    nFreq = length(emData.freqArray)
    freqIdx = collect(Int, 1:nFreq)
    if length(args)>0
        epField = args[1]
        length(args)>1 && ( freqIdx = args[2] )

        isempty(epField) && ( epField = zeros(0) )
        isempty(freqIdx) && ( freqIdx = collect(Int, 1:nFreq) )
    end

    @printf("%-40s", "Forming receiver projection matrix:")
    t = @elapsed rxMat = compEMRxTerm(emData.rxLoc, emMesh)
    @printf("%8.2f %s\n", t, "seconds.")

    # extract things
    nTx   = size(emData.txLoc, 1)
    nRx   = size(emData.rxLoc, 1)
    nFreq = length(freqIdx)

    nx, ny, nz = emMesh.gridSize

    # number of edge variables
    nEx = nx * (ny+1) * (nz+1)
    nEy = (nx+1) * ny * (nz+1)
    nEz = (nx+1) * (ny+1) * nz
    nE  = nEx + nEy + nEz

    # number of face variables
    nFx = (nx+1) * ny * nz
    nFy = nx * (ny+1) * nz
    nFz = nx * ny * (nz+1)
    nF  = nFx + nFy + nFz

    # compute all components of EM fields
    Ex = Array{ComplexF64}(undef, nRx, nTx, nFreq)
    Ey = copy(Ex)
    Ez = copy(Ex)

    # *****
    interpTotalEz = true
    isempty(epField) && (interpTotalEz = false)

    for iFreq = 1:nFreq
        tmp = view(eField, :, :, iFreq)
        Ex[:, :, iFreq] = rxMat.exMat' * tmp
        Ey[:, :, iFreq] = rxMat.eyMat' * tmp

        interpTotalEz && ( tmp += view(epField, :, :, iFreq) )
        Ez[:, :, iFreq] = rxMat.ezMat' * tmp
    end

    if size(bField, 1)>0
        Bx = Array{ComplexF64}(undef, nRx, nTx, nFreq)
        By = copy(Bx)
        Bz = copy(Bx)

        for iFreq = 1:nFreq
            tmpB = view(bField, :, :, iFreq)
            tmpE = view(eField, :, :, iFreq)
            Bx[:, :, iFreq] = rxMat.bxMat' * tmpB + rxMat.bxeMat' * tmpE
            By[:, :, iFreq] = rxMat.byMat' * tmpB + rxMat.byeMat' * tmpE
            Bz[:, :, iFreq] = rxMat.bzMat' * tmpB
        end

    else
        Bx = zeros(0)
        By = zeros(0)
        Bz = zeros(0)
    end


    Exp, Eyp, Ezp, Bxp, Byp, Bzp = getPrimaryFieldsAtRxs(emData, refModel, freqIdx)

    Ex += Exp
    Ey += Eyp

    # if interpTotalEz, Ez is already the total field.
    !interpTotalEz && ( Ez += Ezp )

    if size(bField, 1)>0
        Bx += Bxp
        By += Byp
        Bz += Bzp
    end

    return Ex, Ey, Ez, Bx, By, Bz

end


"""
Interpolate fields to receiver positions using the receiver projection matrix.

Input:
    emMesh    :: EMTensorMesh
    emData    :: MTData
    eField    :: Array   - secondary electric fields at grid sampling points.
    bField    :: Array   - secondary magnetic fields at grid sampling points.
    freqIdx   :: Vector  - Frequency index, for parallel processing (optional).

Output:
    total fields (all six components) at receivers for all freqs, Tx.

"""
function getTotalFieldsAtRxsMT(emMesh::EMTensorMesh, emData::MTData,
         eField::Array{T,3}, bField::Array{T,3}, freqIdx=zeros(0)) where {T<:Complex}

    # deal with all frequencies by default
    if isempty(freqIdx)
        nFreq = length(emData.freqArray)
        freqIdx = collect(Int, 1:nFreq)
    end

    @printf("%-40s", "Forming receiver projection matrix:")
    t = @elapsed rxMat = compEMRxTerm(emData.rxLoc, emMesh)
    @printf("%8.2f %s\n", t, "seconds.")

    # extract things
    nTx   = 2
    nRx   = size(emData.rxLoc, 1)
    nFreq = length(emData.freqArray[freqIdx])

    nx, ny, nz = emMesh.gridSize

    # number of edge variables
    nEx = nx * (ny+1) * (nz+1)
    nEy = (nx+1) * ny * (nz+1)
    nEz = (nx+1) * (ny+1) * nz
    nE  = nEx + nEy + nEz

    # number of face variables
    nFx = (nx+1) * ny * nz
    nFy = nx * (ny+1) * nz
    nFz = nx * ny * (nz+1)
    nF  = nFx + nFy + nFz

    # compute all components of EM fields
    Ex = Array{ComplexF64}(undef, nRx, nTx, nFreq)
    Ey = copy(Ex)
    Ez = zeros(0)

    Bx = Array{ComplexF64}(undef, nRx, nTx, nFreq)
    By = copy(Bx)
    Bz = copy(Bx)


    for iFreq = 1:nFreq
        tmpE = view(eField, :, :, iFreq)
        tmpB = view(bField, :, :, iFreq)

        Ex[:, :, iFreq] = rxMat.exMat' * tmpE        # (nRx, 2)
        Ey[:, :, iFreq] = rxMat.eyMat' * tmpE
        # Ez[:, :, iFreq] = rxMat.ezMat' * tmpE

        Bx[:, :, iFreq] = rxMat.bxMat' * tmpB + rxMat.bxeMat' * tmpE
        By[:, :, iFreq] = rxMat.byMat' * tmpB + rxMat.byeMat' * tmpE
        Bz[:, :, iFreq] = rxMat.bzMat' * tmpB
    end

    return Ex, Ey, Ez, Bx, By, Bz

end


"""
Interpolate fields to receiver positions (an alternative method).
"""
function getTotalFieldsAtRxsMT2(emMesh::EMTensorMesh, emData::MTData,
         eField::Array{T,3}, bField::Array{T,3}, freqIdx=zeros(0)) where {T<:Complex}

    # deal with all frequencies by default
    if isempty(freqIdx)
        nFreq = length(emData.freqArray)
        freqIdx = collect(Int, 1:nFreq)
    end

    # extract things
    rxLoc = emData.rxLoc
    xLen = emMesh.xLen
    yLen = emMesh.yLen
    zLen = emMesh.zLen

    nx, ny, nz = emMesh.gridSize

    xNode = cumsum([0; xLen]) .- emMesh.origin[1]
    yNode = cumsum([0; yLen]) .- emMesh.origin[2]
    zNode = cumsum([0; zLen]) .- emMesh.origin[3]

    # find out on the top of which layer the receivers are located.
    # assume all reveivers have the same z-loc (no topography).
    zRx = rxLoc[1,3]
    zid = findfirst(x -> x<0.01, abs.(zNode .- zRx))

    zLen1 = zLen[zid]

    if emMesh.condType == "isotropy"
        sigma1 = reshape(emMesh.sigma, nx, ny, nz)[:, :, zid]

    elseif emMesh.condType == "anisotropy"
        sigma1 = zeros(Float64, nx, ny, 6)
        sigma1[:, :, 1] = reshape(emMesh.sigma[:,1], nx, ny, nz)[:, :, zid]
        sigma1[:, :, 2] = reshape(emMesh.sigma[:,2], nx, ny, nz)[:, :, zid]
        sigma1[:, :, 3] = reshape(emMesh.sigma[:,3], nx, ny, nz)[:, :, zid]

        if !isempty(emMesh.offsigma)
            sigma1[:, :, 4] = reshape(emMesh.offsigma[:,1], nx, ny, nz)[:, :, zid]
            sigma1[:, :, 5] = reshape(emMesh.offsigma[:,2], nx, ny, nz)[:, :, zid]
            sigma1[:, :, 6] = reshape(emMesh.offsigma[:,3], nx, ny, nz)[:, :, zid]
        end

    end


    nEx = nx * (ny+1) * (nz+1)
    nEy = (nx+1) * ny * (nz+1)
    nEz = (nx+1) * (ny+1) * nz

    nBx = (nx+1) * ny * nz
    nBy = nx * (ny+1) * nz
    nBz = nx * ny * (nz+1)

    nPol  = 2
    nRx   = size(rxLoc, 1)
    nFreq = length(emData.freqArray[freqIdx])

    # compute all components of EM fields
    Ex = Array{ComplexF64}(undef, nRx, nPol, nFreq)
    Ey = copy(Ex)
    Ez = zeros(0)

    Bx = Array{ComplexF64}(undef, nRx, nPol, nFreq)
    By = copy(Bx)
    Bz = copy(Bx)


    for i=1:nFreq

        # 1st interpolate EM fields at grid nodes to receiver locations.
        # nPol = 2
        # EHr = zeros(ComplexF64, nRx, 5, nPol)

        for j=1:nPol

            ex12 = reshape(eField[1:nEx, j, i], nx, ny+1, nz+1)[:, :, zid:zid+1]
            ey12 = reshape(eField[nEx+1:nEx+nEy, j, i], nx+1, ny, nz+1)[:, :, zid:zid+1]
            bz12 = reshape(bField[nBx+nBy+1:end, j, i], nx, ny, nz+1)[:, :, zid:zid+1]
            bx1  = reshape(bField[1:nBx, j, i], nx+1, ny, nz)[:, :, zid]
            by1  = reshape(bField[nBx+1:nBx+nBy, j, i], nx, ny+1, nz)[:, :, zid]

            if emMesh.condType == "isotropy"
                exr, eyr, bxr, byr, bzr = interpFieldsToRx(rxLoc[:, 1:2], xNode, yNode, zLen1, sigma1,
                                              ex12, ey12, bz12, bx1, by1)

            elseif emMesh.condType == "anisotropy"
                ez1 = reshape(eField[nEx+nEy+1:nEx+nEy+nEz, j, i], nx+1, ny+1, nz)[:, :, zid]

                exr, eyr, bxr, byr, bzr = interpFieldsToRxAni(rxLoc[:, 1:2], xNode, yNode, zLen1, sigma1,
                                              ex12, ey12, bz12, bx1, by1, ez1)

            end

            Ex[:, j, i] = exr
            Ey[:, j, i] = eyr
            Bx[:, j, i] = bxr
            By[:, j, i] = byr
            Bz[:, j, i] = bzr

        end # nPol


    end  # nFreq

    return Ex, Ey, Ez, Bx, By, Bz

end



"""
`interpFieldsToRxAni` computes EM fields at receiver locations from fields at grid nodes
for a single frequency with a single polarization mode.

Input:
    rxLoc  :: Array     - receiver locations.
    xNode  :: Vector    - node x-coordinates.
    yNode  :: Vector    - node y-coordinates.
    zLen1  :: Float64   - z-cell size of the receiver layer.
    sigma1 :: Array     - cell conductivity of the receiver layer.
    Ex, Ey, Bz :: Array{ComplexF64, 3}    - grid fields of the receiver layer.
    Ez, Bx, By :: Array{ComplexF64, 2}    - grid fields of the receiver layer.

Output:
    Er, Hr :: Array{ComplexF64,nRx,1}   - Fields at the receiver locations.

"""
function interpFieldsToRxAni(rxLoc::Array{T1,2}, xNode::Vector{T1}, yNode::Vector{T1},
                        zLen1::T1, sigma1::Array{T1,3}, Ex::Array{T2,3},
                        Ey::Array{T2,3}, Bz::Array{T2,3}, Bx::Array{T2,2},
                        By::Array{T2,2}, Ez::Array{T2,2}) where {T1<:Real, T2<:Complex}


    MU0 = 4 * pi * 1e-7

    xLen = diff(xNode)
    yLen = diff(yNode)

    # cell center coordiantes
    xCen = xNode[1:end-1] + xLen / 2
    yCen = yNode[1:end-1] + yLen / 2

    nx = length(xLen)
    ny = length(yLen)

    sig1xx = sigma1[:,:,1]
    sig1yy = sigma1[:,:,2]
    sig1zz = sigma1[:,:,3]
    sig1xy = sigma1[:,:,4]
    sig1xz = sigma1[:,:,5]
    sig1yz = sigma1[:,:,6]

    # First compute Hy at the receiver layer (earth surface or seafloor).
    Hy0 = zeros(ComplexF64, nx, ny+1)

    for ix=1:nx
        Bz0 = vec(Bz[ix, :, 1])
        Bz1 = vec(Bz[ix, :, 2])

        Ex0 = vec(Ex[ix, 2:end-1, 1])
        Ex1 = vec(Ex[ix, 2:end-1, 2])

        Ey0 = 0.5 * ( vec(Ey[ix, :, 1]) + vec(Ey[ix+1, :, 1]) )
        Ey1 = 0.5 * ( vec(Ey[ix, :, 2]) + vec(Ey[ix+1, :, 2]) )

        EzH = 0.5 * ( vec(Ez[ix, 2:end-1]) + vec(Ez[ix+1, 2:end-1]) )

        # quarter Hz (1/4), with length ny
        HzQ = (0.75 * Bz0 + 0.25 * Bz1) / MU0

        # half Hy (1/2), with length ny-1
        # More strictly, an average mu should be used here.
        HyH = vec(By[ix, 2:end-1]) / MU0

        # quarter Ex (1/4), with length ny-1
        ExQ = 0.75 * Ex0 + 0.25 * Ex1

        # quarter Ey collocated with ExQ
        Ey0t = (av(ny-1)*(Ey0./yLen)) .* (yLen[1:end-1].*yLen[2:end]) ./ (av(ny-1)*yLen)
        Ey1t = (av(ny-1)*(Ey1./yLen)) .* (yLen[1:end-1].*yLen[2:end]) ./ (av(ny-1)*yLen)
        EyQ  = 0.75 * Ey0t + 0.25 * Ey1t

        # quarter Ez collocated with ExQ, assume Ez=0 at the earth surface
        EzQ = 0.5 * EzH


        # average conductiviy at vertical edge
        sigtxx = vec(sig1xx[ix,:])
        sigtxy = vec(sig1xy[ix,:])
        sigtxz = vec(sig1xz[ix,:])

        sig1xxV = (av(ny-1)*(sigtxx.*yLen)) ./ (av(ny-1)*yLen)
        sig1xyV = (av(ny-1)*(sigtxy.*yLen)) ./ (av(ny-1)*yLen)
        sig1xzV = (av(ny-1)*(sigtxz.*yLen)) ./ (av(ny-1)*yLen)

        JxQ = sig1xxV.*ExQ + sig1xyV.*EyQ + sig1xzV.*EzQ

        # dHz/dy
        dHzQ = (ddx(ny-1)*HzQ) ./ (av(ny-1)*yLen)

        # Ampre's theorem: dHz/dy - dHy/dz = JxQ
        # where dHy/dz = (HyH-Hy0)/(0.5*zLen1).
        Hy0[ix, 2:end-1] = HyH - (dHzQ - JxQ)*(0.5*zLen1)
        Hy0[ix, 1]   = Hy0[ix, 2]
        Hy0[ix, end] = Hy0[ix, end-1]
    end  # nx


    # Second compute Hx at the receiver layer (earth surface or seafloor).
    Hx0 = zeros(ComplexF64, nx+1, ny)

    for iy=1:ny
        Bz0 = Bz[:, iy, 1]
        Bz1 = Bz[:, iy, 2]

        Ey0 = Ey[2:end-1, iy, 1]
        Ey1 = Ey[2:end-1, iy, 2]

        Ex0 = 0.5 * ( Ex[:, iy, 1] + Ex[:, iy+1, 1] )
        Ex1 = 0.5 * ( Ex[:, iy, 2] + Ex[:, iy+1, 2] )

        EzH = 0.5 * ( Ez[2:end-1, iy] + Ez[2:end-1, iy+1] )

        # quarter Hz (1/4), with length nx
        HzQ = (0.75 * Bz0 + 0.25 * Bz1) / MU0

        # half Hx (1/2), with length nx-1
        # More strictly, an average mu should be used here.
        HxH = Bx[2:end-1, iy] / MU0

        # quarter Ey (1/4), with length nx-1
        EyQ = 0.75 * Ey0 + 0.25 * Ey1

        # quarter Ex collocated with EyQ
        Ex0t = (av(nx-1)*(Ex0./xLen)) .* (xLen[1:end-1].*xLen[2:end]) ./ (av(nx-1)*xLen)
        Ex1t = (av(nx-1)*(Ex1./xLen)) .* (xLen[1:end-1].*xLen[2:end]) ./ (av(nx-1)*xLen)
        ExQ  = 0.75 * Ex0t + 0.25 * Ex1t

        # quarter Ez collocated with EyQ, assume Ez=0 at the earth surface
        EzQ = 0.5 * EzH

        # average conductiviy at vertical edge
        sigtyy = vec(sig1yy[:,iy])
        sigtxy = vec(sig1xy[:,iy])
        sigtyz = vec(sig1yz[:,iy])

        sig1yyV = (av(nx-1)*(sigtyy.*xLen)) ./ (av(nx-1)*xLen)
        sig1xyV = (av(nx-1)*(sigtxy.*xLen)) ./ (av(nx-1)*xLen)
        sig1yzV = (av(nx-1)*(sigtyz.*xLen)) ./ (av(nx-1)*xLen)

        JyQ = sig1xyV.*ExQ + sig1yyV.*EyQ + sig1yzV.*EzQ

        # dHz/dx
        dHzQ = (ddx(nx-1)*HzQ) ./ (av(nx-1)*xLen)

        # Ampre's theorem: dHx/dz - dHz/dx = JyQ,
        # where dHx/dz = (HxH-Hx0)/(0.5*zLen1).
        Hx0[2:end-1, iy] = HxH - (dHzQ + JyQ)*(0.5*zLen1)
        Hx0[1, iy]   = Hx0[2, iy]
        Hx0[end, iy] = Hx0[end-1, iy]
    end  # ny


    # Third interpolate fields to receiver locations (using bilinear interpolation)
    # pre-defined fields at receiver locations
    nRx = size(rxLoc,1)
    Exr = zeros(ComplexF64, nRx)
    Eyr = copy(Exr)
    Bxr = copy(Exr)
    Byr = copy(Exr)
    Bzr = copy(Exr)


    itpMat1 = bilinearInterpMat(rxLoc, xCen, yNode)    # for Ex & Hy
    itpMat2 = bilinearInterpMat(rxLoc, xNode, yCen)    # for Ey & Hx
    itpMat3 = bilinearInterpMat(rxLoc, xCen, yCen)     # for Hz

    Exr = itpMat1' * vec(Ex[:, :, 1])
    Eyr = itpMat2' * vec(Ey[:, :, 1])
    Bxr = itpMat2' * vec(Hx0) * MU0
    Byr = itpMat1' * vec(Hy0) * MU0
    Bzr = itpMat3' * vec(Bz[:, :, 1])

    return Exr, Eyr, Bxr, Byr, Bzr
end


"""
`interpFieldsToRx` is for isotropic case.
"""
function interpFieldsToRx(rxLoc::Array{T1,2}, xNode::Vector{T1}, yNode::Vector{T1},
                        zLen1::T1, sigma1::Array{T1,2}, Ex::Array{T2,3},
                        Ey::Array{T2,3}, Bz::Array{T2,3}, Bx::Array{T2,2},
                        By::Array{T2,2}) where {T1<:Real, T2<:Complex}


    MU0 = 4 * pi * 1e-7

    xLen = diff(xNode)
    yLen = diff(yNode)

    # cell center coordiantes
    xCen = xNode[1:end-1] + xLen / 2
    yCen = yNode[1:end-1] + yLen / 2

    nx = length(xLen)
    ny = length(yLen)

    # First compute Hy at the receiver layer (earth surface or seafloor).
    Hy0 = zeros(ComplexF64, nx, ny+1)

    for iy=1:ny-1
        Bz0L = Bz[:, iy, 1]
        Bz1L = Bz[:, iy, 2]
        Bz0R = Bz[:, iy+1, 1]
        Bz1R = Bz[:, iy+1, 2]

        Ex0 = Ex[:, iy+1, 1]
        Ex1 = Ex[:, iy+1, 2]

        # quarter Hz (1/4)
        HzQL = (0.75 * Bz0L + 0.25 * Bz1L) / MU0
        HzQR = (0.75 * Bz0R + 0.25 * Bz1R) / MU0

        # half Hy (1/2)
        # More strictly, an average mu should be used here.
        HyH = By[:, iy+1] / MU0

        # quarter Ex (1/4)
        ExQ = 0.75 * Ex0 + 0.25 * Ex1

        # average conductiviy at vertical edge
        #sigtL = vec(sigma1[:,iy])
        #sigtR = vec(sigma1[:,iy+1])
        sigma1v = ( sigma1[:,iy]*yLen[iy]+sigma1[:,iy+1]*yLen[iy+1] ) / ( yLen[iy]+yLen[iy+1] )

        # dHz/dy
        dHzQ = 2*(HzQR-HzQL) / ( yLen[iy]+yLen[iy+1] )

        # Ampre's theorem: dHz/dy - dHy/dz = sigma1v.*ExQ,
        # where dHy/dz = (HyH-Hy0)/(0.5*zLen1).
        Hy0[:, iy+1] = HyH - (dHzQ - sigma1v.*ExQ)*(0.5*zLen1)
    end  # ny
    Hy0[:, 1]   = Hy0[:, 2]
    Hy0[:, end] = Hy0[:, end-1]


    # Second compute Hx at the receiver layer (earth surface or seafloor).
    Hx0 = zeros(ComplexF64, nx+1, ny)

    for iy=1:ny
        Bz0 = Bz[:, iy, 1]
        Bz1 = Bz[:, iy, 2]

        Ey0 = Ey[2:end-1, iy, 1]
        Ey1 = Ey[2:end-1, iy, 2]

        # quarter Hz (1/4), with length nx
        HzQ = (0.75 * Bz0 + 0.25 * Bz1) / MU0

        # half Hx (1/2), with length nx-1
        # More strictly, an average mu should be used here.
        HxH = Bx[2:end-1, iy] / MU0

        # quarter Ey (1/4), with length nx-1
        EyQ = 0.75 * Ey0 + 0.25 * Ey1

        # average conductiviy at vertical edge
        sigt = sigma1[:,iy]
        sigma1v = (av(nx-1)*(sigt.*xLen)) ./ (av(nx-1)*xLen)

        # dHz/dx
        dHzQ = (ddx(nx-1)*HzQ) ./ (av(nx-1)*xLen)

        # Ampre's theorem: dHx/dz - dHz/dx = sigma1v.*EyQ,
        # where dHx/dz = (HxH-Hx0)/(0.5*zLen1).
        Hx0[2:end-1, iy] = HxH - (dHzQ + sigma1v.*EyQ)*(0.5*zLen1)
        Hx0[1, iy]   = Hx0[2, iy]
        Hx0[end, iy] = Hx0[end-1, iy]
    end  # ny


    # Third interpolate fields to receiver locations (using bilinear interpolation)
    # pre-defined fields at receiver locations
    nRx = size(rxLoc,1)
    Exr = zeros(ComplexF64, nRx)
    Eyr = copy(Exr)
    Bxr = copy(Exr)
    Byr = copy(Exr)
    Bzr = copy(Exr)


    itpMat1 = bilinearInterpMat(rxLoc, xCen, yNode)    # for Ex & Hy
    itpMat2 = bilinearInterpMat(rxLoc, xNode, yCen)    # for Ey & Hx
    itpMat3 = bilinearInterpMat(rxLoc, xCen, yCen)     # for Hz

    Exr = itpMat1' * vec(Ex[:, :, 1])
    Eyr = itpMat2' * vec(Ey[:, :, 1])
    Bxr = itpMat2' * vec(Hx0) * MU0
    Byr = itpMat1' * vec(Hy0) * MU0
    Bzr = itpMat3' * vec(Bz[:, :, 1])

    return Exr, Eyr, Bxr, Byr, Bzr
end


"""
`getTransFunction` computes MT transfer functions (e.g., impedance, apparent
 resistivity, tipper) from EM fields.

Input:
    freqs           :: Array   - frequencies.
    Ex,Ey,Bx,By,Bz  :: Array   - EM fields at the receiver locations.
    dataType        :: String  - the type of responses (impedance or rho & phase).

Output:
    resp :: Array
"""
function getTransFunction(freqs::Array, Ex::Array{T, 3}, Ey::Array{T, 3},
                          Bx::Array{T, 3}, By::Array{T, 3}, Bz::Array{T, 3},
                          dataType::AbstractString) where {T<:Complex}

    MU0 = 4 * pi * 1e-7

    nFreq = length(freqs)
    nRx   = size(Ex, 1)

    Zxx = Array{ComplexF64}(undef, nRx, 1, nFreq)
    Zxy = copy(Zxx)
    Zyx = copy(Zxx)
    Zyy = copy(Zxx)
    Tzx = copy(Zxx)
    Tzy = copy(Zxx)

    for iFreq = 1:nFreq
        ex1 = Ex[:, 1, iFreq];   ex2 = Ex[:, 2, iFreq]
        ey1 = Ey[:, 1, iFreq];   ey2 = Ey[:, 2, iFreq]
        hx1 = Bx[:, 1, iFreq]/MU0;   hx2 = Bx[:, 2, iFreq]/MU0
        hy1 = By[:, 1, iFreq]/MU0;   hy2 = By[:, 2, iFreq]/MU0
        hz1 = Bz[:, 1, iFreq]/MU0;   hz2 = Bz[:, 2, iFreq]/MU0

        # Impedance:
        #   [Ex1 Ex2     [Zxx Zxy     [Hx1 Hx2
        #    Ey1 Ey2] =   Zyx Zyy] *   Hy1 Hy2]
        # so:
        #   [Zxx Zxy     [Ex1 Ex2     [Hy2 -Hx2
        #    Zyx Zyy] =   Ey1 Ey2] *  -Hy1  Hx1] / det(H)
        detH = hx1 .* hy2 - hy1 .* hx2
        Zxx[:, 1, iFreq] = (ex1 .* hy2 - ex2 .* hy1) ./ detH
        Zxy[:, 1, iFreq] = (ex2 .* hx1 - ex1 .* hx2) ./ detH
        Zyx[:, 1, iFreq] = (ey1 .* hy2 - ey2 .* hy1) ./ detH
        Zyy[:, 1, iFreq] = (ey2 .* hx1 - ey1 .* hx2) ./ detH


        # Tipper:
        #   [Hz1 Hz2] = [Tzx Tzy] * [Hx1 Hx2
        #                            Hy1 Hy2]
        # so:
        #   [Tzx Tzy] = [Hz1 Hz2] * [Hy2 -Hx2
        #                           -Hy1  Hx1] / det(H)
        Tzx[:, 1, iFreq] = (hy2 .* hz1 - hy1 .* hz2) ./ detH
        Tzy[:, 1, iFreq] = (hx1 .* hz2 - hx2 .* hz1) ./ detH
    end


    if dataType == "Impedance"
        resp = hcat(Zxx, Zxy, Zyx, Zyy)

    elseif dataType == "Impedance_Tipper"
        resp = hcat(Zxx, Zxy, Zyx, Zyy, Tzx, Tzy)

    elseif occursin("Rho_Phs", dataType)
        rhoxx = Array{Float64}(undef, nRx, 1, nFreq)
        rhoxy = copy(rhoxx)
        rhoyx = copy(rhoxx)
        rhoyy = copy(rhoxx)
        phsxx = copy(rhoxx)
        phsxy = copy(rhoxx)
        phsyx = copy(rhoxx)
        phsyy = copy(rhoxx)

        for iFreq = 1:nFreq
            omega = 2 * pi * freqs[iFreq]
            Ztmp = Zxx[:, 1, iFreq]
            rhoxx[:, 1, iFreq] = (abs.(Ztmp)).^2 / (omega*MU0)
            phsxx[:, 1, iFreq] = atan.(imag(Ztmp), real(Ztmp)) * 180/pi

            Ztmp = Zxy[:, 1, iFreq]
            rhoxy[:, 1, iFreq] = (abs.(Ztmp)).^2 / (omega*MU0)
            phsxy[:, 1, iFreq] = atan.(imag(Ztmp), real(Ztmp)) * 180/pi

            Ztmp = Zyx[:, 1, iFreq]
            rhoyx[:, 1, iFreq] = (abs.(Ztmp)).^2 / (omega*MU0)
            phsyx[:, 1, iFreq] = atan.(imag(Ztmp), real(Ztmp)) * 180/pi

            Ztmp = Zyy[:, 1, iFreq]
            rhoyy[:, 1, iFreq] = (abs.(Ztmp)).^2 / (omega*MU0)
            phsyy[:, 1, iFreq] = atan.(imag(Ztmp), real(Ztmp)) * 180/pi
        end

        if occursin("Tipper", dataType)
            resp = hcat(rhoxx, phsxx, rhoxy, phsxy, rhoyx, phsyx, rhoyy, phsyy,
                        real(Tzx), imag(Tzx), real(Tzy), imag(Tzy))
        else
            resp = hcat(rhoxx, phsxx, rhoxy, phsxy, rhoyx, phsyx, rhoyy, phsyy)
        end

    end

    return resp
end


#
# `matchData`, it's not a pure forward problem by default, then we need to
# select the data. This one is for CSEM.
#
function matchData(emData::CSEMData, fwdResp::Array{T}) where {T<:Complex}

    # extract things
    txID = emData.txID
    rxID = emData.rxID
    dtID = emData.dtID
	dataType = emData.dataType
    freqID   = emData.freqID

    nTx   = size(emData.txLoc, 1)
    nFreq = length(emData.freqArray)

    nData = length(txID)
    predData = zeros(ComplexF64, nData)
    p = 1
    for iFreq = 1:nFreq

        # predicted at receiver locations
        # first find source indices associated with freq i
        indF  = findall(freqID .== iFreq)

        # check if current frequency exists
        isempty(indF) && continue

        subTxID = txID[indF]
        for iTx = 1:nTx

            # extract source indices
            indTx = findall(subTxID .== iTx)

            # check if current source exists at specified freq
            isempty(indTx) && continue

            subRxID = rxID[indF][indTx]
            subDtID = dtID[indF][indTx]
            nd  = length(subRxID)
            for k = 1:nd
                iRx = subRxID[k]
                iDt = subDtID[k]
                dt  = dataType[iDt]

                if dt == "Ex"
                    predData[p] = fwdResp[iRx, iTx, iFreq]
                elseif dt == "Ey"
                    predData[p] = fwdResp[iRx, nTx + iTx, iFreq]
                elseif dt == "Ez"
                    predData[p] = fwdResp[iRx, 2*nTx + iTx, iFreq]
                elseif dt == "Bx"
                    predData[p] = fwdResp[iRx, 3*nTx + iTx, iFreq]
                elseif dt == "By"
                    predData[p] = fwdResp[iRx, 4*nTx + iTx, iFreq]
                elseif dt == "Bz"
                    predData[p] = fwdResp[iRx, 5*nTx + iTx, iFreq]
                end

                p += 1
            end # k = nd

        end # nTx

    end # nFreq

    return predData

end


#
# `matchData`, it's not a pure forward problem by default, then we need to
# select the data. This one is for MT.
#
function matchData(emData::MTData, fwdResp::Array{T}) where {T<:Union{Real, Complex}}

    # extract things
    rxID     = emData.rxID
    dcID     = emData.dcID
    freqID   = emData.freqID
    dataType = emData.dataType
    dataComp = emData.dataComp

    nFreq = length(emData.freqArray)

    nData = length(rxID)
    if occursin("Impedance", dataType)
        predData = zeros(ComplexF64, nData)
    elseif occursin("Rho_Phs", dataType)
        predData = zeros(Float64, nData)
    end

    p = 1
    for i = 1:nFreq

        indF  = findall(freqID .== i)

        # check if current frequency exists
        isempty(indF) && continue

        subRxID = rxID[indF]
        subDcID = dcID[indF]
        nd  = length(subRxID)

        for j = 1:nd

            iRx = subRxID[j]
            iDt = subDcID[j]
            dt  = dataComp[iDt]

            if occursin("Impedance", dataType)

                if dt == "ZXX"
                    predData[p] = fwdResp[iRx, 1, i]
                elseif dt == "ZXY"
                    predData[p] = fwdResp[iRx, 2, i]
                elseif dt == "ZYX"
                    predData[p] = fwdResp[iRx, 3, i]
                elseif dt == "ZYY"
                    predData[p] = fwdResp[iRx, 4, i]
                end

                if occursin("Tipper", dataType)
                    if dt == "TZX"
                        predData[p] = fwdResp[iRx, 5, i]
                    elseif dt == "TZY"
                        predData[p] = fwdResp[iRx, 6, i]
                    end
                end

            elseif occursin("Rho_Phs", dataType)

                if dt == "RhoXX"
                    predData[p] = fwdResp[iRx, 1, i]
                elseif dt == "PhsXX"
                    predData[p] = fwdResp[iRx, 2, i]
                elseif dt == "RhoXY"
                    predData[p] = fwdResp[iRx, 3, i]
                elseif dt == "PhsXY"
                    predData[p] = fwdResp[iRx, 4, i]
                elseif dt == "RhoYX"
                    predData[p] = fwdResp[iRx, 5, i]
                elseif dt == "PhsYX"
                    predData[p] = fwdResp[iRx, 6, i]
                elseif dt == "RhoYY"
                    predData[p] = fwdResp[iRx, 7, i]
                elseif dt == "PhsYY"
                    predData[p] = fwdResp[iRx, 8, i]
                end

                if occursin("Tipper", dataType)
                    if dt == "RealTZX"
                        predData[p] = fwdResp[iRx, 9, i]
                    elseif dt == "ImagTZX"
                        predData[p] = fwdResp[iRx, 10, i]
                    elseif dt == "RealTZY"
                        predData[p] = fwdResp[iRx, 11, i]
                    elseif dt == "ImagTZY"
                        predData[p] = fwdResp[iRx, 12, i]
                    end
                end

            else
                error("Data type $(dataType) is not supported!")

            end

            p += 1

        end # j = nd

    end # i = nFreq

    return predData

end
