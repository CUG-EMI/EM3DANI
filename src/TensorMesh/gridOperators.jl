#------------------------------------------------------------------------------
# Define routines for discretizing differential operators with respect to the
# tensor mesh.
#------------------------------------------------------------------------------


# Discretize face divergence operator mapping from face center to cell center.
function getFaceDivergence(xLen::Vector, yLen::Vector, zLen::Vector)

    nx = length(xLen)
    ny = length(yLen)
    nz = length(zLen)

    Dx = kron3(spunit(nz), spunit(ny), ddx(nx))
    Dy = kron3(spunit(nz), ddx(ny), spunit(nx))
    Dz = kron3(ddx(nz), spunit(ny), spunit(nx))

    faceDivMat = [Dx Dy Dz]

    faceArea = meshGeoFace(xLen, yLen, zLen)
    volM     = meshGeoVolumeInv(xLen, yLen, zLen)

    faceDivMat  = volM * faceDivMat * faceArea

    return faceDivMat

end


# Discretize nodal gradient operator mapping from node to edge center.
function getNodalGradient(xLen::Vector, yLen::Vector, zLen::Vector)

    nx = length(xLen)
    ny = length(yLen)
    nz = length(zLen)

    Gx = kron3(spunit(nz+1), spunit(ny+1), ddx(nx))
    Gy = kron3(spunit(nz+1), ddx(ny), spunit(nx+1))
    Gz = kron3(ddx(nz), spunit(ny+1), spunit(nx+1))

    gradM = [Gx; Gy; Gz]

    edgeM = meshGeoEdgeInv(xLen, yLen, zLen);

    gradM = edgeM * gradM

    return gradM

end


# Discretize cell gradient operator mapping from cell center to cell face.
function getCellGradient(xLen::Vector, yLen::Vector, zLen::Vector)

    nx = length(xLen)
    ny = length(yLen)
    nz = length(zLen)

    Gx = kron3(spunit(nz), spunit(ny), ddxC2N(nx))
    Gy = kron3(spunit(nz), ddxC2N(ny), spunit(nx))
    Gz = kron3(ddxC2N(nz), spunit(ny), spunit(nx))

    cellGradMat = [Gx; Gy; Gz]

    # get face areas and cell volumes
    faceMat = meshGeoFace(xLen, yLen, zLen)
    volMat  = meshGeoVolume(xLen, yLen, zLen)
    AvCF    = aveCell2Face([nx, ny, nz])
    volMat  = AvCF * diag(volMat)

    disMat   = diag(faceMat) ./ volMat
    cellGradMat = sdiag(disMat) * cellGradMat

    nxFace = size(Gx,1)
    nyFace = size(Gy,1)

    xcGrad = cellGradMat[1:nxFace, :]
    ycGrad = cellGradMat[nxFace+1:nxFace + nyFace, :]
    zcGrad = cellGradMat[nxFace + nyFace + 1 : end, :]

    return cellGradMat, xcGrad, ycGrad, zcGrad

end


# Discretize cell gradient operator with boundary conditions.
function getCellGradientBC(xLen::Vector, yLen::Vector, zLen::Vector,
                           BC::Vector{String})

    nx = length(xLen)
    ny = length(yLen)
    nz = length(zLen)

    #BC = "neumann"
    Dx = kron3(spunit(nz), spunit(ny), ddxCellGradBC(nx,BC[1]))
    Dy = kron3(spunit(nz), ddxCellGradBC(ny,BC[2]), spunit(nx))
    Dz = kron3(ddxCellGradBC(nz,BC[3]), spunit(ny), spunit(nx))

    cellGradMat = [Dx; Dy; Dz]

    # get face areas and cell volumes
    faceMat = meshGeoFace(xLen, yLen, zLen)
    volMat  = meshGeoVolume(xLen, yLen, zLen)
    AvCF    = aveCell2Face([nx, ny, nz])
    volMat  = AvCF * diag(volMat);

    disMat   = diag(faceMat) ./ volMat
    cellGradMat = sdiag(disMat) * cellGradMat

    return cellGradMat

end


# Discretize edge curl operator mapping from cell edge to cell face.
function getEdgeCurl(xLen::Vector, yLen::Vector, zLen::Vector)

    nx = length(xLen)
    ny = length(yLen)
    nz = length(zLen)

    nxFace = (nx+1) * ny * nz
    nyFace = nx * (ny+1) * nz
    nzFace = nx * ny * (nz+1)
    nxEdge = nx * (ny+1) * (nz+1)
    nyEdge = (nx+1) * ny * (nz+1)
    nzEdge = (nx+1) * (ny+1) * nz

    # The curl operator from edges to faces
    Dyz = kron3(-ddx(nz), spunit(ny), spunit(nx+1))
    Dzy = kron3(spunit(nz), -ddx(ny), spunit(nx+1))

    Dxz = kron3(-ddx(nz), spunit(ny+1), spunit(nx))
    Dzx = kron3(spunit(nz), spunit(ny+1), -ddx(nx))

    Dxy = kron3(spunit(nz+1), -ddx(ny), spunit(nx))
    Dyx = kron3(spunit(nz+1), spunit(ny), -ddx(nx))

    ExCurl = [spzeros(nxFace, nxEdge) Dyz  -Dzy]
    EyCurl = [-Dxz spzeros(nyFace, nyEdge)  Dzx]
    EzCurl = [Dxy -Dyx  spzeros(nzFace, nzEdge)]

    edgeCurlM = [ExCurl; EyCurl; EzCurl]

    # face area and edge size
    faceMat = meshGeoFaceInv(xLen, yLen, zLen)
    edgeMat = meshGeoEdge(xLen, yLen, zLen)

    edgeCurlM = faceMat * edgeCurlM * edgeMat

    return edgeCurlM

end


#----------------------------------------------------------
# Averaging operators over mesh
#----------------------------------------------------------

# Compute average matrix mapping from face center to cell center.
function aveFace2Cell(n::Vector{Int})

    AvFCx = kron3(spunit(n[3]), spunit(n[2]), av(n[1]))
    AvFCy = kron3(spunit(n[3]), av(n[2]), spunit(n[1]))
    AvFCz = kron3(av(n[3]), spunit(n[2]), spunit(n[1]))

    aveFCMat = [AvFCx AvFCy AvFCz]

    return aveFCMat

end

# Compute average matrix mapping from cell center to face center.
function aveCell2Face(n::Vector{Int})

    AvCFx = kron3(spunit(n[3]), spunit(n[2]), avcn(n[1]))
    AvCFy = kron3(spunit(n[3]), avcn(n[2]), spunit(n[1]))
    AvCFz = kron3(avcn(n[3]), spunit(n[2]), spunit(n[1]))

    aveCFMat = [AvCFx; AvCFy; AvCFz]

    return aveCFMat

end

# Compute average matrix mapping from cell-center to edge.
function aveCell2Edge(n::Vector{Int})

    Avcex = avcn(n[1])
    Avcey = avcn(n[2])
    Avcez = avcn(n[3])

    AvCEx = kron3(Avcez, Avcey, spunit(n[1]))
    AvCEy = kron3(Avcez, spunit(n[2]), Avcex)
    AvCEz = kron3(spunit(n[3]), Avcey, Avcex)

    aveCEMat = [AvCEx; AvCEy; AvCEz]

    return aveCEMat

end

# Compute average matrix averaging mapping from edge to cell-center.
function aveEdge2Cell(n::Vector{Int})

    Avecx = kron3(av(n[3]), av(n[2]), spunit(n[1]))
    Avecy = kron3(av(n[3]), spunit(n[2]), av(n[1]))
    Avecz = kron3(spunit(n[3]), av(n[2]), av(n[1]))

    aveECMat = [Avecx Avecy Avecz]

    return aveECMat

end

# Compute average matrix mapping from node to cell-center.
function aveNode2Cell(n::Vector{Int})

    Ancx = av(n[1])
    Ancy = av(n[2])
    Ancz = av(n[3])

    aveNCMat = kron3(Ancz, Ancy, Ancx)

    return aveNCMat

end

# Compute average matrix mapping from cell-center to node.
function aveCell2Node(n::Vector{Int})

    Avcnx = avcn(n[1])
    Avcny = avcn(n[2])
    Avcnz = avcn(n[3])

    aveCNMat = kron3(Avcnz, Avcny, Avcnx)

    return aveCNMat
end
