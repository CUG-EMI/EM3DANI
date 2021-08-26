
@testset "TensorMesh" begin

	xLen = [1000.; 800; 600; 200; 200; 500]
	yLen = [1200.; 1000; 500; 200; 100; 100; 100]
	zLen = [2000.; 1000; 800; 500; 100; 50; 100; 500]

	# meshGeometry
	@testset "meshGeometry" begin

		@testset "meshGeoFace" begin
			faceM = meshGeoFace(xLen, yLen, zLen)
			@test size(faceM) == (1154, 1154)
			@test nnz(faceM) == 1154
			@test isapprox(norm(diag(faceM))/1e7, 1.708232127, atol=1e-8)
		end

		@testset "meshGeoFaceInv" begin
			faceMInv = meshGeoFaceInv(xLen, yLen, zLen)
			@test size(faceMInv) == (1154, 1154)
			@test nnz(faceMInv) == 1154
			@test isapprox(norm(diag(faceMInv)), 0.0013713219, atol=1e-8)
		end

		@testset "meshGeoEdge" begin
			edgeM = meshGeoEdge(xLen, yLen, zLen)
			@test size(edgeM) == (1321, 1321)
			@test nnz(edgeM) == 1321
			@test isapprox(norm(diag(edgeM)), 26205.72456, atol=1e-5)
		end

		@testset "meshGeoEdgeInv" begin
			edgeMInv = meshGeoEdgeInv(xLen, yLen, zLen)
			@test size(edgeMInv) == (1321, 1321)
			@test nnz(edgeMInv) == 1321
			@test isapprox(norm(diag(edgeMInv)), 0.2435400377, atol=1e-8)
		end

		@testset "meshGeoVolume" begin
			volM = meshGeoVolume(xLen, yLen, zLen)
			@test size(volM) == (336, 336)
			@test nnz(volM) == 336
			@test isapprox(norm(diag(volM))/1e9, 6.295220806, atol=1e-8)
		end

		@testset "meshGeoVolumeInv" begin
			volMInv = meshGeoVolumeInv(xLen, yLen, zLen)
			@test size(volMInv) == (336, 336)
			@test nnz(volMInv) == 336
			@test isapprox(norm(diag(volMInv))/1e-6, 3.4621206234, atol=1e-8)
		end

	end # meshGeometry


	# gridOperators
	@testset "gridOperators" begin
		gridSize = [length(xLen); length(yLen); length(zLen)]

		@testset "getFaceDivergence" begin
			faceDivMat = getFaceDivergence(xLen, yLen, zLen)
			@test size(faceDivMat) == (336, 1154)
			@test nnz(faceDivMat)  == 2016
		end

		@testset "getNodalGradient" begin
			gradM = getNodalGradient(xLen, yLen, zLen)
			@test size(gradM) == (1321, 504)
			@test nnz(gradM)  == 2642
		end

		@testset "getEdgeCurl" begin
			edgeCurlM = getEdgeCurl(xLen, yLen, zLen)
			@test size(edgeCurlM) == (1154, 1321)
			@test nnz(edgeCurlM)  == 4616
		end

		@testset "aveEdge2Cell" begin
			aveEC = aveEdge2Cell(gridSize)
			@test size(aveEC) == (336, 1321)
			@test nnz(aveEC)  == 4032
		end

		@testset "aveFace2Cell" begin
			aveFC = aveFace2Cell(gridSize)
			@test size(aveFC) == (336, 1154)
			@test nnz(aveFC)  == 2016
		end

		@testset "aveNode2Cell" begin
			aveNC = aveNode2Cell(gridSize)
			@test size(aveNC) == (336, 504)
			@test nnz(aveNC)  == 2688
		end

	end # gridOperators


	# getDiscreteOperators
	@testset "getDiscreteOperators" begin
		emMesh = initEMTensorMesh(xLen, yLen, zLen)
		@test norm(emMesh.xLen - xLen) ≈ 0.
		@test norm(emMesh.yLen - yLen) ≈ 0.
		@test norm(emMesh.zLen - zLen) ≈ 0.
		@test emMesh.gridSize == [6; 7; 8]

		getDiscreteOperators!(emMesh)
		@test size(emMesh.curlM) == (1154, 1321)
		@test nnz(emMesh.curlM)  == 4616
		@test size(emMesh.divM) == (336, 1154)
		@test nnz(emMesh.divM)  == 2016
		@test size(emMesh.gradM) == (1321, 504)
		@test nnz(emMesh.gradM)  == 2642
		@test size(emMesh.aveEC) == (336, 1321)
		@test nnz(emMesh.aveEC)  == 4032
		@test size(emMesh.aveFC) == (336, 1154)
		@test nnz(emMesh.aveFC)  == 2016
		@test size(emMesh.aveNC) == (336, 504)
		@test nnz(emMesh.aveNC)  == 2688
		@test size(emMesh.volM) == (336, 336)
		@test nnz(emMesh.volM) == 336
	end # getDiscreteOperators



end
