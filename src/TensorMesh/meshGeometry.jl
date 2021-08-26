#------------------------------------------------------------------------------
# Define routines for computing tensor mesh geometries
#------------------------------------------------------------------------------

# Compute face area matrix of the 3D grid.
function meshGeoFace(xLen::Vector, yLen::Vector, zLen::Vector)

    nx = length(xLen)
    ny = length(yLen)
    nz = length(zLen)

    f1 = kron3( sdiag(zLen), sdiag(yLen), spunit(nx+1) )
    f2 = kron3( sdiag(zLen), spunit(ny+1), sdiag(xLen) )
    f3 = kron3( spunit(nz+1), sdiag(yLen), sdiag(xLen) )

    faceM = blockdiag(blockdiag(f1, f2), f3)

    return faceM
end


# Compute inverse of face area matrix of the 3D grid.
function meshGeoFaceInv(xLen::Vector, yLen::Vector, zLen::Vector)

    nx = length(xLen)
    ny = length(yLen)
    nz = length(zLen)

    f1 = kron3( sdiag(1 ./ zLen), sdiag(1 ./ yLen), spunit(nx+1) )
    f2 = kron3( sdiag(1 ./ zLen), spunit(ny+1), sdiag(1 ./ xLen) )
    f3 = kron3( spunit(nz+1), sdiag(1 ./ yLen), sdiag(1 ./ xLen) )

    faceMInv = blockdiag(blockdiag(f1, f2), f3)

    return faceMInv

end


# Compute edge length matrix of the 3D grid.
function meshGeoEdge(xLen::Vector, yLen::Vector, zLen::Vector)

    nx = length(xLen)
    ny = length(yLen)
    nz = length(zLen)

    e1 = kron3( spunit(nz+1), spunit(ny+1), sdiag(xLen) )
    e2 = kron3( spunit(nz+1), sdiag(yLen), spunit(nx+1) )
    e3 = kron3( sdiag(zLen), spunit(ny+1), spunit(nx+1) )

    # convert to sparse matrix
    edgeM = blockdiag(blockdiag(e1, e2), e3)

    return edgeM
end


# Compute inverse of edge length matrix of the 3D grid.
function meshGeoEdgeInv(xLen::Vector, yLen::Vector, zLen::Vector)

    nx = length(xLen)
    ny = length(yLen)
    nz = length(zLen)

    e1 = kron3( spunit(nz+1), spunit(ny+1), sdiag(1 ./ xLen) )
    e2 = kron3( spunit(nz+1), sdiag(1 ./ yLen), spunit(nx+1) )
    e3 = kron3( sdiag(1 ./ zLen), spunit(ny+1), spunit(nx+1) )

    # convert to sparse matrix
    edgeMInv = blockdiag(blockdiag(e1, e2), e3)

    return edgeMInv

end


# Compute volume matrix of the 3D grid.
function meshGeoVolume(xLen::Vector, yLen::Vector, zLen::Vector)

    volM = kron3(sdiag(zLen), sdiag(yLen), sdiag(xLen))

    return volM
end


# Compute inverse of volume matrix of the 3D grid.
function meshGeoVolumeInv(xLen::Vector, yLen::Vector, zLen::Vector)

    volMInv = kron3(sdiag(1 ./ zLen), sdiag(1 ./ yLen), sdiag(1 ./ xLen))

    return volMInv
end
