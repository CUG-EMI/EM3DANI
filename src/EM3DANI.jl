#=
######## ##     ##  #######  ########     ###    ##    ## ####
##       ###   ### ##     ## ##     ##   ## ##   ###   ##  ##
##       #### ####        ## ##     ##  ##   ##  ####  ##  ##
######   ## ### ##  #######  ##     ## ##     ## ## ## ##  ##
##       ##     ##        ## ##     ## ######### ##  ####  ##
##       ##     ## ##     ## ##     ## ##     ## ##   ###  ##
######## ##     ##  #######  ########  ##     ## ##    ## ####

-------------------------------------------------------------------------------
 A package for isotropic/anisotropic 3D forward modeling of frequency-domain
 electromagnetic (CSEM and MT) data.

 - by Ronghua Peng and Bo Han, China University of Geosciences, Wuhan, 2019.
-------------------------------------------------------------------------------
=#

__precompile__()

module EM3DANI

include("EMUtils/EMUtils.jl")
include("TensorMesh/TensorMesh.jl")
include("EMFileIO/EMFileIO.jl")
include("EM1DUtils/EM1DUtils.jl")
include("LinearSolver/LinearSolver.jl")
include("EMFwdSolver/EMFwdSolver.jl")


end
