#
# Solve the linear equation using direct solvers. Only MUMPS and MKL PARDISO
# have been implemented so far.
#
function directSolver(A::SparseMatrixCSC, rhs::AbstractArray, parm::DirectSolverParm)

    # extract things
    solverName = parm.solverName
    nThread    = parm.nThread
    sym        = parm.sym
    ooc        = parm.ooc
    saveFac    = parm.saveFac

    #---> MUMPS
    if solverName ==:mumps
        @printf("%-40s", "   - MUMPS factorization phase:")
        t = @elapsed Ainv = factorMUMPS(A, sym, ooc)
        @printf("%8.2f %s\n", t, "seconds.")

        n    = size(rhs, 1)
        nRhs = size(rhs, 2)
        x  = zeros(ComplexF64, n, nRhs)

        if maximum(abs.(rhs)) == 0.0
            println("All source terms are zeros.")
        else
            @printf("%-40s", "   - MUMPS solution phase:")
            t = @elapsed x = applyMUMPS(Ainv, rhs)
            @printf("%8.2f %s\n", t, "seconds.")
        end

        if !saveFac
            destroyMUMPS(Ainv)
            Ainv = []
        end

    #---> MKL Pardiso
    elseif solverName ==:mklpardiso

        Ainv = MKLPardisoSolver()

        set_msglvl!(Ainv, Pardiso.MESSAGE_LEVEL_OFF)
        set_nprocs!(Ainv, nThread)

        # you can solve it in a single call
        #@time x = solve(Ainv, A, rhs, :N)

        # it would be better to do it step by step
        set_matrixtype!(Ainv, Pardiso.COMPLEX_SYM)
        pardisoinit(Ainv)

        # out-of-core
        if ooc>0
            set_iparm!(Ainv, 3, 1)
            set_iparm!(Ainv, 60, 2)
        end

        A_pardiso = get_matrix(Ainv, A, :N)

        # Analyze the matrix and compute a symbolic factorization &
        # compute the numeric factorization.
        @printf("%-40s", "   - MKLPardiso factorization phase:")
        set_phase!(Ainv, Pardiso.ANALYSIS_NUM_FACT)
        t = @elapsed pardiso(Ainv, A_pardiso, rhs)
        @printf("%8.2f %s\n", t, "seconds.")

        # Compute the solutions X using the symbolic factorization.
        @printf("%-40s", "   - MKLPardiso solution phase:")
        set_phase!(Ainv, Pardiso.SOLVE_ITERATIVE_REFINE)
        x = similar(rhs)
        t = @elapsed pardiso(Ainv, x, A_pardiso, rhs)
        @printf("%8.2f %s\n", t, "seconds.")

        if !saveFac
            # Free the PARDISO data structures.
            set_phase!(Ainv, Pardiso.RELEASE_ALL)
            pardiso(Ainv)
            Ainv = []
        end


    else
        error("$solverName is not supported!")
    end


    return x, Ainv

end
