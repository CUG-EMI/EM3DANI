################################################################################
# Defines some utility routines for sparse matrix manipulation, 1D finite-
# difference approximation, linear interpolation, parallel computing, etc.
################################################################################
module EMUtils

using  SparseArrays, Distributed
using  LinearAlgebra


include("interpUtils.jl")
include("sparseUtils.jl")
include("parallelUtils.jl")
include("eulerRotUtils.jl")

end
