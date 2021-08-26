################################################################################
# Define routines to solve CSEM/MT 1D isotropic/anisotropic forward modeling
# problem. It contains interface to the fortran code Dipole1D.
################################################################################
module EM1DUtils

using  SparseArrays, Distributed
using  LinearAlgebra


include("dipole1D.jl")
include("mt1DFwd.jl")


end # module
