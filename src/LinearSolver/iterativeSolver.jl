#
# Solve the linear equation using iterative solvers. BiCGStab and QMR
# have been implemented so far.
#
function iterativeSolver(Ke::SparseMatrixCSC, rhs::AbstractArray,
                    parm::IterativeSolverParm, iterLogID::IOStream, PC::Function)


    # extract things
    iterMethod = parm.iterMethod
    maxIter    = parm.maxIter
    tol        = parm.tol
    out        = parm.out

    # (n, nRhs) = size(rhs)
    n    = size(rhs, 1)
    nRhs = size(rhs, 2)
    xt   = zeros(ComplexF64, n, nRhs)

    Af(x) = Ke*x

    if iterMethod == :bicgstb
        #Af(x) = Ke*x
        for j=1:nRhs
            iterTime = @elapsed begin
                @printf(iterLogID, "%s %3d\n","Pol No.:", j)
                @printf(iterLogID, "%4s\t%7s\n","iter","relres")

                xt[:,j], flag, resLast, iter, resvec = bicgstb(Af, rhs[:,j], tol=tol,
                                                       maxIter=maxIter, M1=PC, out=out)

                for k=2:length(resvec)
                    @printf(iterLogID, "%4d %11.2e\n", k-1, resvec[k])
                end

            end  # @elapsed

            @printf(iterLogID, "%s %g %s\n","elapsed time:", iterTime, "seconds.")
        end

    elseif iterMethod == :qmr
        #Af(x) = Ke*x
        for j=1:nRhs
            iterTime = @elapsed begin
                @printf(iterLogID, "%s %3d\n","Pol No.:", j)
                @printf(iterLogID, "%4s\t%7s\n","iter","relres")

                xt[:,j], flag, resLast, iter, resvec = qmr(Af, rhs[:,j], tol=tol,
                                                       maxIter=maxIter, M=PC, out=out)

                for k=2:length(resvec)
                    @printf(iterLogID, "%4d %11.2e\n", k-1, resvec[k])
                end

            end  # @elapsed

            @printf(iterLogID, "%s %g %s\n","elapsed time:", iterTime, "seconds.")
        end

    end # iterMethod == ?

    return xt

end

# include("qmr.jl")
