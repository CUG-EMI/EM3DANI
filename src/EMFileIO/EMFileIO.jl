################################################################################
# Define routines for data and model input and output.
################################################################################
module EMFileIO

using  Printf, SparseArrays
using  EM3DANI.EMUtils
using  EM3DANI.TensorMesh

export EMData, CSEMData, MTData
#-------------------------------------------------------------------------------

abstract type EMData end


# Struct `CSEMData` and `MTData` contain data configurations, e.g. receiver
# locations, frequencies, data types, but does not contain data (and error) itself.
mutable struct CSEMData{T<:Real}  <:  EMData

    txLoc     :: Array{T}                   # source location array
    rxLoc     :: Array{T}                   # receiver location array
    dpLen     :: T                          # dipole length
    phaseCon  :: String                     # phase convention
    freqArray :: Vector{T}                  # frequency array
    dataType  :: Vector{String}             # data type
    txID      :: Vector{Int}                # source index
    rxID      :: Vector{Int}                # receiver index
    freqID    :: Vector{Int}                # frequency index
    dtID      :: Vector{Int}                # datatype index

end


mutable struct MTData{T<:Real}  <:  EMData

    rxLoc     :: Array{T}                   # receiver location array
    phaseCon  :: String                     # phase convention
    freqArray :: Vector{T}                  # frequency array
    dataType  :: String                     # data type
    dataComp  :: Vector{String}             # data component
    rxID      :: Vector{Int}                # receiver index
    freqID    :: Vector{Int}                # frequency index
    dcID      :: Vector{Int}                # datacomp index

end


include("readEMData.jl")
include("writeEMData.jl")
include("readEMModel.jl")
include("writeEMFwdResp.jl")
include("writeEdgeFields.jl")
include("readEdgeFields.jl")

end
