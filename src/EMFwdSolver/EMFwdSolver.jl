################################################################################
# Define the top-level functions for solving 3D EM forward problem, including:
# forming the coefficient matrix and RHS, incorporating boundary conditions,
# solving the linear system, interpolating fields to receiver positions and
# computing the required form of forward responses.
################################################################################
module EMFwdSolver

using Printf, SparseArrays
using Distributed
using LinearAlgebra
using EM3DANI.EM1DUtils
using EM3DANI.EMFileIO
using EM3DANI.TensorMesh
using EM3DANI.EMUtils
using EM3DANI.LinearSolver
# using ProgressMeter

using MUMPS: MUMPSfactorization
using Pardiso: MKLPardisoSolver


export RxProjMat, RefModel

#-------------------------------------------------------------------------------
# struct `RxProjMat` defines interpolation matrix for receivers.
mutable struct RxProjMat

    # electric field
    exMat::SparseMatrixCSC                  # x component
    eyMat::SparseMatrixCSC                  # y component
    ezMat::SparseMatrixCSC                  # z component

    # magnetic field
    bxMat::SparseMatrixCSC                  # x component
    byMat::SparseMatrixCSC                  # y component
    bzMat::SparseMatrixCSC                  # z component

    # additional
    bxeMat::SparseMatrixCSC
    byeMat::SparseMatrixCSC

end # struct


# struct `RefModel` encapsulates reference model.
mutable struct RefModel{T<:Real}

    sig1D::Vector{T}
    depth1D::Vector{T}

end # mutable struct


include("solveEM3DFwd.jl")
include("parsolveEM3DFwd.jl")
include("compPrimaryField.jl")
include("getAniMassMatrix.jl")
include("interpField.jl")
include("dataFunctional.jl")
include("boundaryFields.jl")

end # EMFwdSolver
