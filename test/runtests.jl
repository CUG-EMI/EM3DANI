#push!(LOAD_PATH, "/home/username/code")     # absolute path, for example
push!(LOAD_PATH, pwd() * "/../..")        # relative path

using LinearAlgebra, SparseArrays
using Distributed, Printf
using Test
using EM3DANI.EM1DUtils
#using EM3DANI.EMFileIO
using EM3DANI.TensorMesh
import EM3DANI.TensorMesh: meshGeoFace, meshGeoFaceInv, meshGeoEdge, meshGeoEdgeInv,
                           meshGeoVolume, meshGeoVolumeInv, getFaceDivergence,
                           getNodalGradient, getEdgeCurl, aveEdge2Cell,
                           aveFace2Cell, aveNode2Cell
using EM3DANI.EMUtils
using EM3DANI.LinearSolver


#using MUMPS: MUMPSfactorization
#using Pardiso: MKLPardisoSolver


@testset "EM3DANI" begin

include("testEMUtils.jl")
include("testTensorMesh.jl")
include("testLinearSolver.jl")
include("testEM1DUtils.jl")

end
