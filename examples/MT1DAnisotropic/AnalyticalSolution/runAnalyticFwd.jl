#------------------------------------------------------------------------------
# script to perform MT 1D anisotropic forward modeling by using the EM1DUtils
# module.
#------------------------------------------------------------------------------
#push!(LOAD_PATH, "/home/username/code")     # absolute path, for example
push!(LOAD_PATH, pwd() * "/../../../..")     # relative path
using EM3DANI.EMUtils: eulerRotation!
using EM3DANI.EM1DUtils
using Printf
#------------------------------------------------------------------------------

include("readEM1DAniParam.jl")
include("writeMT1DAniResp.jl")

paramfile = "Pek1DModel_31Freqs.in"

printstyled("Reading 1D anisotropic parameter file $(paramfile) ...\n", color=:blue)
@time freqs, sigma, zNode = readEM1DAniParam(paramfile)

Z, Rho, Phs  = mt1DAnalyticImpAniPek(freqs, sigma, zNode)

respfile = "Pek1DModel_31Freqs.resp"
writeMT1DAniResp(respfile, Z, Rho, Phs)

println("=== Finishing forward problem ===")
