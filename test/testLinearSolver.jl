
@testset "LinearSolver" begin

	# generate a complex symmetric matrix
	A = sprandn(100,100,.1) + 1im*(sprandn(100,100,.1))
	A  = triu(A) + transpose(triu(A,1)) + SparseMatrixCSC((10.0+10im)I, 100, 100)
	rhs = complex(randn(100))


	# iterative solver
	@testset "iterative solver" begin
		lsParm = IterativeSolverParm()

		D  = Vector(diag(A))
		PC = x -> x./D
		lsParm.tol = 1e-7
		lsParm.maxIter = 50

		@testset "qmr" begin
			lsParm.iterMethod = :qmr
			lsParm.out = 0
			iterLogID = open("iterQMR.log", "w")
			x1 = iterativeSolver(A, rhs, lsParm, iterLogID, PC)
			@test norm(A*x1[:, 1]-rhs)/norm(rhs) < lsParm.tol
		end

		@testset "bicgstb" begin
			lsParm.iterMethod = :bicgstb
			lsParm.out = 0
			iterLogID = open("iterBiCGSTB.log", "w")
			x2 = iterativeSolver(A, rhs, lsParm, iterLogID, PC)
			@test norm(A*x2[:, 1]-rhs)/norm(rhs) < lsParm.tol
		end

	end # iterative solver


	# direct solver
	@testset "direct solver" begin
		lsParm = DirectSolverParm()
		lsParm.ooc = 1                      # 0=in-core, 1=out-of-core
		tol = 1e-15

		@testset "mumps" begin
			lsParm.solverName = :mumps
			x3, = directSolver(A, rhs, lsParm)
			@test norm(A*x3[:, 1]-rhs)/norm(rhs) < tol
		end

		@testset "mklpardiso" begin
			lsParm.solverName = :mklpardiso
			x4, = directSolver(A, rhs, lsParm)
			@test norm(A*x4[:, 1]-rhs)/norm(rhs) < tol
		end

	end # direct solver

end
