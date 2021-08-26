
@testset "EMUtils" begin

	# interpUtils
	@testset "interpUtils" begin

		@testset "linearInterp" begin
			x = [1.0, 2, 3, 4, 5, 6, 7]

			p1 = -5.
			indL1, indR1, wL1, wR1 = linearInterp(p1, x)
			@test (indL1, indR1) == (1, 1); @test wL1 ≈ 0.5; @test wR1 ≈ 0.5

			p2 = 1.
			indL2, indR2, wL2, wR2 = linearInterp(p2, x)
			@test (indL2, indR2) == (1, 1); @test wL2 ≈ 0.5; @test wR2 ≈ 0.5

			p3 = 3.2
			indL3, indR3, wL3, wR3 = linearInterp(p3, x)
			@test (indL3, indR3) == (3, 4); @test wL3 ≈ 0.8; @test wR3 ≈ 0.2

			p4 = 4.
			indL4, indR4, wL4, wR4 = linearInterp(p4, x)
			@test (indL4, indR4) == (3, 4); @test wL4 ≈ 0.0; @test wR4 ≈ 1.0

			p5 = 7.
			indL5, indR5, wL5, wR5 = linearInterp(p5, x)
			@test (indL5, indR5) == (6, 7); @test wL5 ≈ 0.0; @test wR5 ≈ 1.0

			p6 = 9.
			indL6, indR6, wL6, wR6 = linearInterp(p6, x)
			@test (indL6, indR6) == (7, 7); @test wL6 ≈ 0.5; @test wR6 ≈ 0.5

		end

		@testset "linearInterpMat" begin
			x = [1.0, 2, 3, 4, 5, 6, 7]
			points = [-5., 1., 3.2, 4., 7., 9.]
			inds, weights = linearInterpMat(points, x, false)
			@test inds[1] == (1, 1)
			@test inds[3] == (3, 4)
			@test inds[5] == (6, 7)
			interpMat = linearInterpMat(points, x)
			@test size(interpMat) == (7, 6)
			@test nnz(interpMat) == 7
			#nzIdx = findall(!iszero, interpMat)
			#@test nzIdx[1] == CartesianIndex(1, 1)
		end

		@testset "bilinearInterpMat" begin
			x = collect(-10.0:10.0)
			y = collect(-12.0:2.0:12.0)
			points = [-5. 1.; 3.2  4.; 7. 9.]
			inds, weights = bilinearInterpMat(points, x, y, false)
			@test inds[1][1] == (5, 7)
			@test inds[2][1] == (14, 8)
			@test inds[3][1] == (17, 11)
			interpMat = bilinearInterpMat(points, x, y)
			@test size(interpMat) == (273, 3)
			@test nnz(interpMat) == 6
		end

		@testset "trilinearInterpMat" begin
			x = collect(-10.0:10.0)
			y = collect(-12.0:2.0:12.0)
			z = collect(-5.3:2.0:15.0)
			points = [-5. 1. 3.2;  4. 7. 9.]
			inds, weights = trilinearInterpMat(points, x, y, z, false)
			@test inds[1][1] == (5, 7, 5)
			@test inds[2][1] == (14, 10, 8)
			interpMat = trilinearInterpMat(points, x, y, z)
			@test size(interpMat) == (3003, 2)
			@test nnz(interpMat) == 8
		end

	end # interpUtils

	# sparseUtils
	@testset "sparseUtils" begin
		A = spunit(3)
		@test nnz(A) == 3
		@test Array(diag(A)) == [1.0; 1.0; 1.0]

		B = spdiags([1.0; 2; 3])
		@test nnz(B) == 3
		@test Array(diag(B)) == [1.0; 2; 3]

		C = sdiag([1.0; 2; 3])
		@test nnz(C) == 3
		@test Array(diag(C)) == [1.0; 2; 3]

		D = ddx(3)
		@test nnz(D) == 6
		@test Array(D) == [-1.0 1 0 0; 0 -1 1 0; 0 0 -1 1]

		E = ddxC2N(3)
		@test nnz(E) == 6
		@test Array(E) == [1.0 0 0; -1 1 0; 0 -1 1; 0 0 1]

		F = av(3)
		@test nnz(F) == 6
		@test Array(F) == [0.5 0.5 0 0; 0 0.5 0.5 0; 0 0 0.5 0.5]

		G = avcn(3)
		@test nnz(G) == 6
		@test Array(G) == [1.0 0 0; 0.5 0.5 0; 0 0.5 0.5; 0 0 1]

		@testset "meshgrid" begin
			x = collect(1.:3);  y = collect(1.:5);
			X, Y = meshgrid(x, y)
			@test X == [1. 2 3; 1 2 3; 1 2 3; 1 2 3; 1 2 3]
			@test Y == [1. 1 1; 2 2 2; 3 3 3; 4 4 4; 5 5 5]
			z = collect(7.:8)
			X, Y, Z = meshgrid(x, y, z)
			@test Z[:, :, 1] == [7. 7 7; 7 7 7; 7 7 7; 7 7 7; 7 7 7]
			@test Z[:, :, 2] == [8. 8 8; 8 8 8; 8 8 8; 8 8 8; 8 8 8]
		end

	end # sparseUtils

	# eulerRotUtils
	@testset "eulerRotUtils" begin
		sigma    = [1.  1 10;  1 10 1; 1.  1 10;  1 10 1]
		aniAngle = [0. 90  0; 90  0 0; 0. 30  0; 30  0 0]
		sigma, offsigma = eulerRotation!(sigma, aniAngle)
		@test sigma[1:2, :] == [1. 10 1; 10 1 1]
		@test sigma[3, 1] == 1.
		@test sigma[4, 3] == 1.
		@test offsigma[1:2, :] == [0. 0 0; 0 0 0]
		@test offsigma[3, 1] ≈ 0.
		@test offsigma[3, 2] ≈ 0.
		@test !(offsigma[3, 3] ≈ 0.)
		@test !(offsigma[4, 1] ≈ 0.)
		@test offsigma[4, 2] ≈ 0.
		@test offsigma[4, 3] ≈ 0.
	end # eulerRotUtils

	# parallelUtils
	@testset "parallelUtils" begin

		@testset "assignTasks" begin
			idx1 = assignTasks(8, 3)
			@test collect(idx1[1]) == [1; 2; 3]
			@test collect(idx1[2]) == [4; 5; 6]
			@test collect(idx1[3]) == [7; 8]
			idx2 = assignTasks(8, 3, 2)
			@test collect(idx2[1]) == [1; 4; 7]
			@test collect(idx2[2]) == [2; 5; 8]
			@test collect(idx2[3]) == [3; 6]
		end

		@testset "initChannel" begin
			a = [1; 2; 3]
			chanRef = initChannel(a)
			@test isready(chanRef)
			@test take!(chanRef) == a
			@test !isready(chanRef)

			# foo = ()->[1; 2; 3]
			# f(x, y) = x + y
			#chanRef2 = initChannel(f, (2,3))
			#@test isready(chanRef2)
			#@test take!(chanRef2) == 5
			#@test !isready(chanRef2)
		end

		@testset "initRemoteChannel" begin
			#addprocs(1)
			a = [1; 2; 3]
			rRef = initRemoteChannel(a)
			@test isready(rRef[1])
			@test take!(rRef[1]) == a
			@test !isready(rRef[1])

		end


	end

end
