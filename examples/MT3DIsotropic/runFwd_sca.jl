#------------------------------------------------------------------------------
# script to perform 3D EM isotropic/anisotropic forward modeling using a
# scattered/secondary-field approach
#------------------------------------------------------------------------------

# Adding the parent directory of the EM3DANI package directory to LOAD_PATH,
# you can either use an absolute or a relative path:
#push!(LOAD_PATH, "/home/username/code")     # absolute path, for example
push!(LOAD_PATH, pwd() * "/../../..")        # relative path

using Printf
using EM3DANI.EMFileIO
using EM3DANI.EMFwdSolver
using EM3DANI.LinearSolver
#------------------------------------------------------------------------------
# begin timing
totaltime = @elapsed begin

println("Getting started ...")

# Setting the number of threads (only useful for the direct solver)
ENV["OMP_NUM_THREADS"] = 8
ENV["MKL_NUM_THREADS"] = 8

# Problem type, "csem" or "mt"
probType = "mt"


#----------------------- read data file and model file ------------------------#
# data file
datafile = "commemi3d2_RhoPhs.dat"
@printf("%-60s", "Reading data file $(datafile):")
t = @elapsed emData = readEMData(datafile, probType)
@printf("%8.4f %s\n", t, "seconds.")

# model file
modfile = "commemi3d2.mod"
@printf("%-60s", "Reading model file $(modfile):")
t = @elapsed emMesh = readEMModel(modfile, sigAir=1e-8)
@printf("%8.4f %s\n", t, "seconds.")

# Since we are using the secondary-field approach, a reference model need
# to be given
refModel = RefModel(zeros(2), zeros(2))
refModel.sig1D   = [1e-8, 1e-2]
refModel.depth1D = [-emMesh.origin[3], 0.]

# Parameter setting for linear solvers
# direct solvers
lsParm = DirectSolverParm()
#lsParm.solverName = :mklpardiso

# iterative solvers
#lsParm = IterativeSolverParm()
#lsParm.iterMethod = :qmr
#lsParm.tol = 1e-7
#lsParm.prec = :aphi


# Calling the (sequential) forward modeling function
fwdResp, = solveEM3DFwd(emMesh, emData, refModel, lsParm)


# Output the result
respfile = "commemi3d2_RhoPhs_sca.resp"
println("Writing out the forward responses ...")
writeEMFwdResp(respfile, emData, fwdResp)


println("=== Finishing forward problem ===")

end  # @elapsed
@printf("%s %10.2f %s\n", "========== Total elapsed time:", totaltime, "seconds. ==========")
