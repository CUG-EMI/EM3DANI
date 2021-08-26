#------------------------------------------------------------------------------
# script to perform 3D EM isotropic/anisotropic forward modeling using a
# total-field approach
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
datafile = "31Freqs_RhoPhs.dat"
@printf("%-60s", "Reading data file $(datafile):")
t = @elapsed emData = readEMData(datafile, probType)
@printf("%8.4f %s\n", t, "seconds.")

# model file
modfile = "Pek1DModel.mod"
@printf("%-60s", "Reading model file $(modfile):")
t = @elapsed emMesh = readEMModel(modfile, sigAir=1e-8)
@printf("%8.4f %s\n", t, "seconds.")

# Although we are using the total-field approach, the inputting argument
# "refModel" need to be constructed, and its value will be neglected.
refModel = RefModel(zeros(0), zeros(0))

# Parameter setting for linear solvers
# direct solvers
lsParm = DirectSolverParm()
lsParm.solverName = :mumps
#lsParm.solverName = :mklpardiso

# iterative solvers
#lsParm = IterativeSolverParm()
#lsParm.iterMethod = :qmr
#lsParm.tol = 1e-7
#lsParm.prec = :aphi


# Calling the (sequential) forward modeling function
# Note that the last argument "secondary" must be set to be false for the
# total-field approach
fwdResp, = solveEM3DFwd(emMesh, emData, refModel, lsParm, secondary=false)

# Output the result
respfile = "31Freqs_RhoPhs_Total_MUMPS.resp"
println("Writing out the forward responses ...")
writeEMFwdResp(respfile, emData, fwdResp)


println("=== Finishing forward problem. ===")

end  # @elapsed
@printf("%s %10.2f %s\n", "========== Total elapsed time:", totaltime, "seconds. ==========")
