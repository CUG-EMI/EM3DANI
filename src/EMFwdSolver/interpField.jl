# export trilinearInterpEH, trilinearInterpJzAni
# export shantInterpJz, shantInterpBx, shantInterpBy
# export getRxSigmaZ, getEdgeSigma, locateDipole

"""
`trilinearInterpEH` gets electric fields at receiver locations using
 trilinear interpolation scheme.

 Input:
    compFlag: component flag, 1,2,3,4,5,6 for Ex, Ey, Ez, Bx, By, Bz, respectively.
"""
function trilinearInterpEH(compFlag::Int, rxLoc::Array{T,2}, emMesh::EMTensorMesh,
                           inJz=true) where {T<:Real}

    # grid node coordiantes
    xNode = cumsum([0; emMesh.xLen]) .- emMesh.origin[1]
    yNode = cumsum([0; emMesh.yLen]) .- emMesh.origin[2]
    zNode = cumsum([0; emMesh.zLen]) .- emMesh.origin[3]

    # cell center coordiantes
    xCen = xNode[1:end-1] + emMesh.xLen / 2
    yCen = yNode[1:end-1] + emMesh.yLen / 2
    zCen = zNode[1:end-1] + emMesh.zLen / 2

    # number of grid nodes
    nx = length(xNode)
    ny = length(yNode)
    nz = length(zNode)

    nr = size(rxLoc, 1)

    # for E-fields and B-fields, respectively
    if compFlag >=1 && compFlag <=3
        # number of edges
        nEx = (nx-1) * ny * nz
        nEy = nx * (ny-1) * nz
        nEz = nx * ny * (nz-1)
        nE  = nEx + nEy + nEz

        spmEx = spzeros(nEx, nr)
        spmEy = spzeros(nEy, nr)
        spmEz = spzeros(nEz, nr)
        rxMat = spzeros(nE, nr)

    elseif compFlag >=4 && compFlag <=6
        # number of faces
        nBx = nx * (ny-1) * (nz-1)
        nBy = (nx-1) * ny * (nz-1)
        nBz = (nx-1) * (ny-1) * nz
        nB  = nBx + nBy + nBz

        spmBx = spzeros(nBx, nr)
        spmBy = spzeros(nBy, nr)
        spmBz = spzeros(nBz, nr)
        rxMat = spzeros(nB, nr)

    else
       error("Undifined E-H component flag: $compFlag !")
    end


    if compFlag == 1
        tmpMat = trilinearInterpMat(rxLoc, xCen, yNode, zNode)
        rxMat  = vcat(tmpMat, spmEy, spmEz)

    elseif compFlag == 2
        tmpMat = trilinearInterpMat(rxLoc, xNode, yCen, zNode)
        rxMat  = vcat(spmEx, tmpMat, spmEz)

    elseif compFlag == 3
        tmpMat = trilinearInterpMat(rxLoc, xNode, yNode, zCen)
        rxMat  = vcat(spmEx, spmEy, tmpMat)

        # interpolate Jz rather than Ez
        if inJz
            edgeSigma = getEdgeSigma(emMesh)

            if size(emMesh.sigma, 2) > 1    # anisotropic
                sigzz = emMesh.sigma[:, 3]
            else
                sigzz = copy(emMesh.sigma)
            end

            rxSigmaZ  = getRxSigmaZ(sigzz, rxLoc, xNode, yNode, zNode)
            rxMat = ( sparse(Diagonal(1 ./ rxSigmaZ)) * rxMat' * sparse(Diagonal(edgeSigma)) )'
            rxMat = copy(rxMat)
        end

    elseif compFlag == 4
        tmpMat = trilinearInterpMat(rxLoc, xNode, yCen, zCen)
        rxMat  = vcat(tmpMat, spmBy, spmBz)

    elseif compFlag == 5
        tmpMat = trilinearInterpMat(rxLoc, xCen, yNode, zCen)
        rxMat  = vcat(spmBx, tmpMat, spmBz)

    elseif compFlag == 6
        tmpMat = trilinearInterpMat(rxLoc, xCen, yCen, zNode)
        rxMat  = vcat(spmBx, spmBy, tmpMat)
    end

    return rxMat

end


"""
`shantInterpJz` interpolates the normal electric fields at receiver locations
  using the scheme proposed by Shantsev et al.(GJI,2015).

Input:
    rxLoc
    emMesh

"""
function shantInterpJz(rxLoc::Array{T,2}, emMesh::EMTensorMesh) where {T<:Real}

    # make sure that rxLoc has only a single row
    # if size(rxLoc,1) != 1
    #    error("rxLoc has more than one row in funtion interpJz!")
    # end

    # grid location
    xNode = cumsum([0; emMesh.xLen]) .- emMesh.origin[1]
    yNode = cumsum([0; emMesh.yLen]) .- emMesh.origin[2]
    zNode = cumsum([0; emMesh.zLen]) .- emMesh.origin[3]

    # mesh information
    nx = length(xNode)
    ny = length(yNode)
    nz = length(zNode)

    if size(emMesh.sigma, 2) > 1    # anisotropic
        sigzz = emMesh.sigma[:, 3]
    else
        sigzz = copy(emMesh.sigma)
    end
    sigma = reshape(sigzz, nx-1, ny-1, nz-1)

    # get number of edges
    nEx = (nx-1) * ny * nz
    nEy = nx * (ny-1) * nz
    nEz = nx * ny * (nz-1)
    nE  = nEx + nEy + nEz

    nr = size(rxLoc, 1)
    rxMat = spzeros(nE, nr)

    zCen = zNode[1:end-1] + diff(zNode) / 2

    inds2D, weights2D = bilinearInterpMat(rxLoc, xNode, yNode, false)
    inds3D,           = trilinearInterpMat(rxLoc, xNode, yNode, zCen, false)


    for i=1:nr
        indxL = inds3D[i][1][1]
        indyL = inds3D[i][1][2]
        indzL = inds3D[i][1][3]

        # get conductivity ratio
        sigma1 = sigma[indxL,indyL,indzL]
        sigma2 = sigma[indxL,indyL,indzL+1]
        sigmaR = sigma2 / sigma1

        # get z- cordinate of the interface
        zIF = zNode[indzL+1]
        z0  = zCen[indzL]
        z1  = zCen[indzL+1]
        zp  = rxLoc[i,3]

        sz = sigmaR * (z1-zIF)
        wzL = (sz + zIF-zp) / (sz + zIF - z0)
        wzR = (zp-z0) / (sz + zIF - z0)

        w1 = weights2D[i,1] * wzL
        w2 = weights2D[i,2] * wzL
        w3 = weights2D[i,3] * wzL
        w4 = weights2D[i,4] * wzL
        w5 = weights2D[i,1] * wzR
        w6 = weights2D[i,2] * wzR
        w7 = weights2D[i,3] * wzR
        w8 = weights2D[i,4] * wzR

        weights = [w1; w2; w3; w4; w5; w6; w7; w8]

        #
        idx = zeros(Int, 8)
        nz  = length(zCen)

        for j = 1:8
            idx[j] = nx * ny * (inds3D[i][j][3]-1) + nx * (inds3D[i][j][2]-1) + inds3D[i][j][1]
            idx[j] += nEx + nEy
        end

        rxMat[:, i] = sparsevec(idx, weights, nE)
    end

    edgeSigma = getEdgeSigma(emMesh)
    rxSigmaZ  = getRxSigmaZ(sigzz, rxLoc, xNode, yNode, zNode)
    rxMat = ( sparse(Diagonal(1 ./ rxSigmaZ)) * rxMat' * sparse(Diagonal(edgeSigma)) )'
    rxMat = copy(rxMat)

    return rxMat

end



"""
`shantInterpBx` interpolates the tangential magnetic field Bx at receiver locations
  using the scheme proposed by Shantsev et al.(GJI,2015).

Input:
    rxLoc
    emMesh

"""
function shantInterpBx(rxLoc::Array{T,2}, emMesh::EMTensorMesh) where {T<:Real}

    # make sure that rxLoc has only a single row
    # if size(rxLoc,1) != 1
    #    error("rxLoc has more than one row in funtion shantInterpBx!")
    # end

    # grid location
    xNode = cumsum([0; emMesh.xLen]) .- emMesh.origin[1]
    yNode = cumsum([0; emMesh.yLen]) .- emMesh.origin[2]
    zNode = cumsum([0; emMesh.zLen]) .- emMesh.origin[3]

    # mesh information
    nx = length(xNode)
    ny = length(yNode)
    nz = length(zNode)

    if size(emMesh.sigma, 2) > 1    # anisotropic
        sigyy = emMesh.sigma[:, 2]
    else
        sigyy = copy(emMesh.sigma)
    end
    sigma = reshape(sigyy, nx-1, ny-1, nz-1)

    # get number of edges
    nEx = (nx-1) * ny * nz
    nEy = nx * (ny-1) * nz
    nEz = nx * ny * (nz-1)
    nE  = nEx + nEy + nEz

    nr = size(rxLoc, 1)
    BxEMat = spzeros(nE, nr)

    yCen = yNode[1:end-1] + diff(yNode) / 2
    zCen = zNode[1:end-1] + diff(zNode) / 2

    inds2D, weights2D = bilinearInterpMat(rxLoc, xNode, yCen, false)
    inds3D,           = trilinearInterpMat(rxLoc, xNode, yCen, zCen, false)


    for i=1:nr
        indxL = inds3D[i][1][1]
        indyL = inds3D[i][1][2]
        indzL = inds3D[i][1][3]

        # get conductivity difference
        sigma1 = sigma[indxL,indyL,indzL]
        sigma2 = sigma[indxL,indyL,indzL+1]
        sigmaD = sigma2 - sigma1

        # get z- cordinate of the interface
        zIF = zNode[indzL+1]
        z0  = zCen[indzL]
        z1  = zCen[indzL+1]
        zp  = rxLoc[i,3]

        delta = -(z1-zIF)*(zp-z0)/(z1-z0)*sigmaD

        weights = vec( weights2D[i,:] * delta )


        # for Ey
        idx = zeros(Int, 4)

        zidex = inds3D[i][5][3]
        na = nx * (ny-1) * (zidex-1)

        for j = 1:4
            idx[j] = nx * (inds2D[i][j][2]-1) + inds2D[i][j][1]
            idx[j] += (nEx + na)
        end

        BxEMat[:, i] = sparsevec(idx, weights, nE)
    end

    return BxEMat

end


"""
`shantInterpBy` interpolates the tangential magnetic field By at receiver locations
  using the scheme proposed by Shantsev et al.(GJI,2015).

Input:
    rxLoc
    emMesh

"""
function shantInterpBy(rxLoc::Array{T,2}, emMesh::EMTensorMesh) where {T<:Real}

    # make sure that rxLoc has only a single row
    # if size(rxLoc,1) != 1
    #    error("rxLoc has more than one row in funtion shantInterpBy!")
    # end

    # grid location
    xNode = cumsum([0; emMesh.xLen]) .- emMesh.origin[1]
    yNode = cumsum([0; emMesh.yLen]) .- emMesh.origin[2]
    zNode = cumsum([0; emMesh.zLen]) .- emMesh.origin[3]

    # mesh information
    nx = length(xNode)
    ny = length(yNode)
    nz = length(zNode)

    if size(emMesh.sigma, 2) > 1    # anisotropic
        sigxx = emMesh.sigma[:, 1]
    else
        sigxx = copy(emMesh.sigma)
    end
    sigma = reshape(sigxx, nx-1, ny-1, nz-1)

    # get number of edges
    nEx = (nx-1) * ny * nz
    nEy = nx * (ny-1) * nz
    nEz = nx * ny * (nz-1)
    nE  = nEx + nEy + nEz

    nr = size(rxLoc, 1)
    ByEMat = spzeros(nE, nr)

    xCen = xNode[1:end-1] + diff(xNode) / 2
    zCen = zNode[1:end-1] + diff(zNode) / 2

    inds2D, weights2D = bilinearInterpMat(rxLoc, xCen, yNode, false)
    inds3D,           = trilinearInterpMat(rxLoc, xCen, yNode, zCen, false)


    for i=1:nr
        indxL = inds3D[i][1][1]
        indyL = inds3D[i][1][2]
        indzL = inds3D[i][1][3]

        # get conductivity difference
        sigma1 = sigma[indxL,indyL,indzL]
        sigma2 = sigma[indxL,indyL,indzL+1]
        sigmaD = sigma2 - sigma1

        # get z- cordinate of the interface
        zIF = zNode[indzL+1]
        z0  = zCen[indzL]
        z1  = zCen[indzL+1]
        zp  = rxLoc[i,3]

        delta = (z1-zIF)*(zp-z0)/(z1-z0)*sigmaD

        weights = vec( weights2D[i,:] * delta )


        # for Ey
        idx = zeros(Int, 4)

        zidex = inds3D[i][5][3]
        na = (nx-1) * ny * (zidex-1)

        for j = 1:4
            idx[j] = (nx-1) * (inds2D[i][j][2]-1) + inds2D[i][j][1]
            idx[j] += na
        end

        ByEMat[:, i] = sparsevec(idx, weights, nE)
    end

    return ByEMat

end


"""
`trilinearInterpJz` gets electric fields at receiver locations using
 trilinear interpolation scheme. General anisotropy is taken into account.

"""
function trilinearInterpJzAni(rxLoc::Array{Float64,2}, emMesh::EMTensorMesh)

    #
    AveEC    = emMesh.aveEC
    Vol      = emMesh.volM
    gridSize = emMesh.gridSize

    sigxx = emMesh.sigma[:, 1]
    sigyy = emMesh.sigma[:, 2]
    sigzz = emMesh.sigma[:, 3]

    nx, ny, nz = gridSize

    # grid node coordiantes
    xNode = cumsum([0; emMesh.xLen]) .- emMesh.origin[1]
    yNode = cumsum([0; emMesh.yLen]) .- emMesh.origin[2]
    zNode = cumsum([0; emMesh.zLen]) .- emMesh.origin[3]

    # cell center coordiantes
    xCen = xNode[1:end-1] + emMesh.xLen / 2
    yCen = yNode[1:end-1] + emMesh.yLen / 2
    zCen = zNode[1:end-1] + emMesh.zLen / 2

    nEx = nx*(ny+1)*(nz+1)
    nEy = (nx+1)*ny*(nz+1)
    nEz = (nx+1)*(ny+1)*nz
    nE  = nEx + nEy + nEz

    nr = size(rxLoc, 1)

    spmEx = spzeros(nEx, nr)
    spmEy = spzeros(nEy, nr)
    spmEz = spzeros(nEz, nr)
    rxMat = spzeros(nE, nr)


    tmpMat = trilinearInterpMat(rxLoc, xNode, yNode, zCen)
    rxMat  = vcat(spmEx, spmEy, tmpMat)

    # interpolate Jz rather than Ez
    AveX = AveEC[:,1:nEx]
    AveY = AveEC[:,nEx+1:nEx+nEy]
    AveZ = AveEC[:,nEx+nEy+1:end]

    Msig = vcat(AveX'*(Vol*(sigxx)), AveY'*(Vol*(sigyy)), AveZ'*(Vol*(sigzz)))
    Msig = sparse(Diagonal(Msig))

    if !isempty(emMesh.offsigma)
        sigxy = emMesh.offsigma[:, 1]
        sigxz = emMesh.offsigma[:, 2]
        sigyz = emMesh.offsigma[:, 3]
        MsigAni = getOffDiagMassMatrix(gridSize, Vol, sigxy, sigxz, sigyz)
        Msig += MsigAni
    end

    edgeVol = AveEC' * Vector( diag(Vol) )
    Msig = sparse(Diagonal(1 ./ edgeVol)) * Msig

    rxSigmaZ  = getRxSigmaZ(sigzz, rxLoc, xNode, yNode, zNode)
    rxMat = ( sparse(Diagonal(1 ./ rxSigmaZ)) * rxMat' * Msig )'

    return rxMat

end


#
# `getRxSigmaZ` gets the conductivities along z-direction at receiver sites.
#
function getRxSigmaZ(sigma::Vector{T}, rxLoc::Array{T,2}, xNode::Vector{T},
                     yNode::Vector{T}, zNode::Vector{T}) where {T<:Real}

    nr = size(rxLoc, 1)
    rxSigmaZ = zeros(eltype(sigma), nr)

    nx = length(xNode)-1;  ny = length(yNode)-1;  nz = length(zNode)-1

    for i = 1:nr
        indx = locateSegment1D(rxLoc[i, 1], xNode)
        indy = locateSegment1D(rxLoc[i, 2], yNode)
        indz = locateSegment1D(rxLoc[i, 3], zNode)
        ind = (indz-1)*nx*ny + (indy-1)*nx + indx
        rxSigmaZ[i] = sigma[ind]
    end

    return rxSigmaZ
end


#
# `getEdgeSigma(emMesh)` abtains edge conductivities.
#
function getEdgeSigma(emMesh::EMTensorMesh)

    #
    AveEC    = emMesh.aveEC
    Vol      = emMesh.volM
    gridSize = emMesh.gridSize

    if size(emMesh.sigma, 2) > 1    # anisotropic
        sigxx = emMesh.sigma[:, 1]
        sigyy = emMesh.sigma[:, 2]
        sigzz = emMesh.sigma[:, 3]
    else
        sigxx = copy(emMesh.sigma)
        sigyy = copy(emMesh.sigma)
        sigzz = copy(emMesh.sigma)
    end

    nx, ny, nz = gridSize
    nEx = nx*(ny+1)*(nz+1)
    nEy = (nx+1)*ny*(nz+1)
    nEz = (nx+1)*(ny+1)*nz

    #
    AveX = AveEC[:,1:nEx]
    AveY = AveEC[:,nEx+1:nEx+nEy]
    AveZ = AveEC[:,nEx+nEy+1:end]

    edgeSigma = vcat(AveX'*(Vol*(sigxx)), AveY'*(Vol*(sigyy)), AveZ'*(Vol*(sigzz)))
    edgeVol   = AveEC' * Vector( diag(Vol) )

    edgeSigma = edgeSigma ./ edgeVol

    return edgeSigma

end


#
# locateDipole finds the location of dipole to the nearest cell center of the grid
#
function locateDipole(point::Array{T,2}, xLoc::Vector{T}, yLoc::Vector{T},
                      zLoc::Vector{T}) where {T<:Real}

    location = ones(Int,3)

    # x location
    val,ind = findmin(abs(point[1] .- xLoc))
    location[1] = ind

    # y location
    val,ind = findmin(abs(point[2] .- yLoc))
    location[2] = ind

    # z location
    val,ind = findmin(abs(point[3] .- zLoc))
    location[3] = ind

    return location

end
