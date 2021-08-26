################################################################################
# Define routines related to the Yee's tensor product grid.
################################################################################
module TensorMesh

using  SparseArrays, Distributed
using  LinearAlgebra, Printf
using  EM3DANI.EMUtils

export EMTensorMesh, initEMTensorMesh, getDiscreteOperators!

#------------------------------------------------------------------------------
mutable struct EMTensorMesh{T<:Real}

    xLen::Vector{T}                    # width of cells in x direction
    yLen::Vector{T}                    # width of cells in y direction
    zLen::Vector{T}                    # width of cells in z direction
    airLayer::Vector{T}                # air layer
    seaLayer::Vector{T}                # sea layer
    gridSize::Vector{Int}              # cell numbers in x,y,z direction
    origin::Vector{T}                  # origin of mesh

    # conductivity type: isotropy, anisotropy
    condType::AbstractString
    sigma::Array{T}
    offsigma::Array{T}

    volM::SparseMatrixCSC            # Volume of cells
    curlM::SparseMatrixCSC           # edge curl
    gradM::SparseMatrixCSC           # nodal gradient
    divM::SparseMatrixCSC            # face divergence
    aveEC::SparseMatrixCSC             # averaging mapping from edge to cell center
    aveFC::SparseMatrixCSC             # averaging mapping from face to cell center
    aveNC::SparseMatrixCSC             # averaging mapping from node to cell center

end


# Initialize EMTensorMesh
function initEMTensorMesh(xLen::Vector, yLen::Vector, zLen::Vector, x0=zeros(3))

    n = [length(xLen); length(yLen); length(zLen)]
    eMat  = zeros(0)
    epMat = spzeros(0, 0)

    mesh = EMTensorMesh(xLen, yLen, zLen, eMat, eMat, n, x0, "", eMat, eMat,
                        epMat, epMat, epMat, epMat, epMat, epMat, epMat)

    return mesh

end


# Set up the discrete differential operators
function getDiscreteOperators!(emMesh::EMTensorMesh)
    xLen = emMesh.xLen
    yLen = emMesh.yLen
    zLen = emMesh.zLen
    gridSize = emMesh.gridSize

    emMesh.curlM  = getEdgeCurl(xLen, yLen, zLen)
    emMesh.divM   = getFaceDivergence(xLen, yLen, zLen)
    emMesh.gradM  = getNodalGradient(xLen, yLen, zLen)
    emMesh.aveEC  = aveEdge2Cell(gridSize)
    emMesh.aveFC  = aveFace2Cell(gridSize)
    emMesh.aveNC  = aveNode2Cell(gridSize)
    emMesh.volM   = meshGeoVolume(xLen, yLen, zLen)

    return emMesh
end

include("meshGeometry.jl")
include("gridOperators.jl")

end
