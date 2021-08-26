################################################################################
# Define solver wrapper for solving linear systems.
# Currently available solvers include:
# MUMPS, MKL PARDISO (direct solvers),
# QMR, BiCGStab (iterative solvers).
################################################################################
module LinearSolver

using SparseArrays, Distributed
using LinearAlgebra, Printf

#------------------------------------------------------------------------------
using MUMPS
using Pardiso
using KrylovMethods

#-------------------------------------------------------------------------------
export LinearSolverParm
export DirectSolverParm, IterativeSolverParm
export directSolver, iterativeSolver

abstract type LinearSolverParm end

# parameters for direct solvers.
mutable struct DirectSolverParm  <:  LinearSolverParm
    solverName::Symbol        # :mumps, :mklpardiso.
    nThread::Int              # number of thread, only used by mklpardiso
    sym::Int                  # 0=unsymmetric, 1=symm. pos def, 2=general symmetric
    ooc::Int                  # 0=in-core, 1=out-of-core
    saveFac::Bool
    DirectSolverParm() = new(:mumps, 1, 1, 0, false)   # Constructor
end


# parameters for iterative solvers.
mutable struct IterativeSolverParm  <:  LinearSolverParm
    iterMethod::Symbol        # iterative method to use
    prec::Symbol              # type of preconditioning (:aphi,:ssor,...)
    # precM::Function           # preconditioner function handle (M \ x)
    maxIter::Int              # maximum number of iterations
    tol::Real                 # error tolerance
    out::Int                  # flag for output
    # Constructor
    IterativeSolverParm() = new(:bicgstb, :aphi, 1000, 1e-6, 2)
end


include("directSolver.jl")
include("iterativeSolver.jl")

end # module
