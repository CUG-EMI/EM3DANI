
@testset "EM1DUtils" begin

	# mt1DFwd
	@testset "mt1DFwd" begin

		# Isotropic case
		@testset "mt1DAnalyticFields" begin
			sigma = 1 ./ [100.; 1000; 1]
			zNode = [0.; 2000; 6000; 10000]
			freqs = [0.001; 0.01; 0.1; 1]
			z = Array{ComplexF64}(undef, length(freqs), 1)
			j = 0
			for j = 1:length(freqs)
				e, h = mt1DAnalyticFields(freqs[j], sigma, zNode, compH=true)
				z[j] = e[1] / h[1]
			end

			z0 = [6.29710529e-05; 2.05275773e-04; 1.043833786e-03; 2.18509224349e-02] + 1im *
			     [1.09909433e-04; 6.67093284e-04; 5.203915598e-03; 3.04114729295e-02]

			@test norm(z-z0)/norm(z0) < 1e-9
		end


		# Anisotropic case (Pek model)
		@testset "mt1DAnalyticAni" begin
			sigma    = 1.0 ./ [1e4 1e4 1e4; 200. 2e4  200.; 1e3 2e3 1e3; 100 100 100]
			aniAngle = [0. 0  0; 15 0 0; -75. 0  0; 0  0 0]
			sigma, offsigma = eulerRotation!(sigma, aniAngle)
			sigma = hcat(sigma, offsigma)
			zNode = 1e3 * [0.; 10; 28; 128; 150]

			freqs = [0.0001; 0.01; 1; 100]

			# mt1DAnalyticImpAniPek
			Z1, (), = mt1DAnalyticImpAniPek(freqs, sigma, zNode[1:end-1])

			# mt1DAnalyticImpAni
			Z2, (), = mt1DAnalyticImpAni(freqs, sigma, zNode[1:end-1])

			# cross-check
			for j = 1:length(freqs)
				@test norm( Z1[:, j] - Z2[:, j] ) / norm( Z1[:, j] ) < 1e-8
			end

			# mt1DAnalyticFieldsAni
			eTop1 = [1.0+0im; 0]     # Ex-Hy
			eTop2 = [0; 1.0+0im]     # Ey-Hx
			for j = 1:length(freqs)
				ex1, ey1, ez1, hx1, hy1 = mt1DAnalyticFieldsAni(freqs[j], sigma, zNode, eTop1, compH=true)
				ex2, ey2, ez2, hx2, hy2 = mt1DAnalyticFieldsAni(freqs[j], sigma, zNode, eTop2, compH=true)

				detH = hx1 .* hy2 - hy1 .* hx2
				Zxx = (ex1 .* hy2 - ex2 .* hy1) ./ detH
				Zxy = (ex2 .* hx1 - ex1 .* hx2) ./ detH
				Zyx = (ey1 .* hy2 - ey2 .* hy1) ./ detH
				Zyy = (ey2 .* hx1 - ey1 .* hx2) ./ detH

				@test isapprox(Zxy[1], Z1[2, j])
				@test isapprox(Zyx[1], Z1[3, j])
			end

		end

	end # mt1DFwd


	# dipole1D
	@testset "dipole1D" begin

		# canonical 1D marine model
		sig1D   = [1e-8, 3.3, 1.0]
		depth1D = [-1e4, 0, 1000.]

		# survey geometry
		freqs = [0.25]
		txLoc = [0. 0 950  90 0; 0. 0 950  0 0]  # inline and broadside
		rxLoc = [0. 2000  1000; 0. 5000  1000; 0. 8000  1000]

		@testset "point dipole" begin
			ex1, ey1, ez1 = comp1DCSEM(txLoc, rxLoc, freqs, sig1D, depth1D)

			ex2 = comp1DCSEM(txLoc, rxLoc, freqs, sig1D, depth1D, 1)
			ey2 = comp1DCSEM(txLoc, rxLoc, freqs, sig1D, depth1D, 2)
			ez2 = comp1DCSEM(txLoc, rxLoc, freqs, sig1D, depth1D, 3)

			@test isapprox(norm(ex1[:, 2, 1]), 5.7171e-12, atol=1e-15)
			@test isapprox(norm(ex2[:, 2, 1]), 5.7171e-12, atol=1e-15)
			@test norm(ex1[:, 2, 1] - ex2[:, 2, 1]) / norm(ex1[:, 2, 1]) < 1e-5

			@test isapprox(norm(ey1[:, 1, 1]), 2.5607e-12, atol=1e-15)
			@test isapprox(norm(ey2[:, 1, 1]), 2.5607e-12, atol=1e-15)
			@test norm(ey1[:, 1, 1] - ey2[:, 1, 1]) / norm(ey1[:, 1, 1]) < 1e-4

			@test isapprox(norm(ez1[:, 1, 1]), 1.4239e-12, atol=1e-15)
			@test isapprox(norm(ez2[:, 1, 1]), 1.4239e-12, atol=1e-15)
			@test norm(ez1[:, 1, 1] - ez2[:, 1, 1]) / norm(ez1[:, 1, 1]) < 1e-4
		end


		@testset "finite-length dipole" begin
			dipLen = 20.0
			nIntPts = 15
			exL = comp1DCSEM(txLoc, rxLoc, dipLen, nIntPts, freqs, sig1D, depth1D, 1)
			eyL = comp1DCSEM(txLoc, rxLoc, dipLen, nIntPts, freqs, sig1D, depth1D, 2)
			ezL = comp1DCSEM(txLoc, rxLoc, dipLen, nIntPts, freqs, sig1D, depth1D, 3)

			@test isapprox(norm(exL[:, 2, 1]), 5.7169e-12, atol=1e-15)
			@test isapprox(norm(eyL[:, 1, 1]), 2.56096e-12, atol=1e-15)
			@test isapprox(norm(ezL[:, 1, 1]), 1.42398e-12, atol=1e-15)
		end

	end # dipole1D

end
