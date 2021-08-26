export readEMModel

"""
Read model parameters from file.

Input:
    modfile   :: String     - the name of the model file.
    sigAir    :: Real       - air conductivity.
    sigWater  :: Real       - seawater conductivity.

Output:
    emMesh    :: EMTensorMesh - model parameters and mesh properties.

"""
function readEMModel(modfile::String; sigAir=1e-8, sigWater=3.3)

    if isfile(modfile)
        fid = open(modfile, "r")
    else
        error("$(modfile) does not exist, please try again.")
    end

    # Initialization
    emp = zeros(0)
    nx   = [];  ny   = [];  nz   = []
    xLen = [];  yLen = [];  zLen = []
    nAir = 0;   airLayer = emp
    nSea = 0;   seaLayer = emp
    origin  = []

    # default form of model parameter: isotropic, linear, conductivity
    condType = "isotropy"
    resType  = "conductivity"
    modType  = "linear"

    sigma  = emp
    sigmax = emp;  sigmay = emp;  sigmaz = emp
    strike = emp;  dip = emp;     slant  = emp
    nBlock = []


    while !eof(fid)

        cline = strip(readline(fid))

        # ignore empty lines and comments (preceded with #)
        isempty(cline)  && continue
        cline[1] == '#' && continue

        # block at x-axis
        if occursin("NX", cline)
            tmp = split(cline)
            nx  = parse(Int, tmp[end])
            nd  = 0
            xLen = zeros(nx)

            while nd < nx
                cline = strip(readline(fid))
                cline = split(cline)
                num = length(cline)
                for i = 1:num
                    nd = nd + 1
                    xLen[nd] = parse(Float64, cline[i])
                end
            end

        # block at y-axis
        elseif occursin("NY", cline)
            tmp = split(cline)
            ny  = parse(Int, tmp[end])
            nd  = 0
            yLen = zeros(ny)

            while nd < ny
                cline = strip(readline(fid))
                cline = split(cline)
                num = length(cline)
                for i = 1:num
                    nd = nd + 1
                    yLen[nd] = parse(Float64, cline[i])
                end
            end

        # block at z-axis
        elseif occursin("NZ", cline)
            tmp = split(cline)
            nz  = parse(Int, tmp[end])
            nd  = 0
            zLen = zeros(nz)

            while nd < nz
                cline = strip(readline(fid))
                cline = split(cline)
                num = length(cline)
                for i = 1:num
                    nd = nd + 1
                    zLen[nd] = parse(Float64, cline[i])
                end
            end

            nBlock = nx * ny * nz

        # air layer
        elseif occursin("NAIR", cline)
            tmp = split(cline)
            nAir  = parse(Int, tmp[end])
            nd  = 0
            airLayer = zeros(nAir)

            while nd < nAir
                cline = strip(readline(fid))
                cline = split(cline)
                num = length(cline)
                for i = 1:num
                    nd = nd + 1
                    airLayer[nd] = parse(Float64, cline[i])
                end
            end

        elseif occursin("NSEA", cline)
            tmp = split(cline)
            nSea = parse(Int, tmp[end])
            nd  = 0
            seaLayer = zeros(nSea)

            while nd < nSea
                cline = strip(readline(fid))
                cline = split(cline)
                num = length(cline)
                for i = 1:num
                    nd = nd + 1
                    seaLayer[nd] = parse(Float64, cline[i])
                end
            end

        elseif occursin("Anisotropy Type", cline)
            tmp = split(cline)
            condType = lowercase(tmp[end])

        # resistivity type: resistivity or conductivity
        elseif occursin("Resistivity Type", cline)
            tmp = split(cline)
            resType = lowercase(tmp[end])

        # model type: linear or logorithmic
        elseif occursin("Model Type", cline)
            tmp = split(cline)
            modType = lowercase(tmp[end])


        # origin
        elseif occursin("Origin", cline)
            tmp = split(cline)
            origin = zeros(3)
            origin[1] = parse(Float64, tmp[end-2])
            origin[2] = parse(Float64, tmp[end-1])
            origin[3] = parse(Float64, tmp[end])

        elseif condType == "isotropy"
            if occursin("sigma:", cline)
                sigma = zeros(nBlock)
                nd = 0
                while nd < nBlock
                    cline = strip(readline(fid))
                    cline = split(cline)
                    num = length(cline)
                    for i = 1:num
                        nd = nd + 1
                        sigma[nd] = parse(Float64, cline[i])
                    end
                end
            end

        elseif condType == "anisotropy"

            if occursin("sigmax", cline)
                sigmax = zeros(nBlock)
                nd = 0
                while nd < nBlock
                    cline = strip(readline(fid))
                    cline = split(cline)
                    num = length(cline)
                    for i = 1:num
                        nd = nd + 1
                        sigmax[nd] = parse(Float64,cline[i])
                    end
                end

            elseif occursin("sigmay", cline)
                sigmay = zeros(nBlock)
                nd = 0
                while nd < nBlock
                    cline = strip(readline(fid))
                    cline = split(cline)
                    num = length(cline)
                    for i = 1:num
                        nd = nd + 1
                        sigmay[nd] = parse(Float64,cline[i])
                    end
                end

            elseif occursin("sigmaz", cline)
                sigmaz = zeros(nBlock)
                nd = 0
                while nd < nBlock
                    cline = strip(readline(fid))
                    cline = split(cline)
                    num = length(cline)
                    for i = 1:num
                        nd = nd + 1
                        sigmaz[nd] = parse(Float64,cline[i])
                    end
                end # while

            # strike angle
            elseif occursin("strike", cline)
                strike = zeros(nBlock)
                nd = 0
                while nd < nBlock
                    cline = strip(readline(fid))
                    cline = split(cline)
                    num = length(cline)
                    for i = 1:num
                        nd = nd + 1
                        strike[nd] = parse(Float64,cline[i])
                    end
                end # while

            # dip angle
            elseif occursin("dip", cline)
                dip = zeros(nBlock)
                nd  = 0
                while nd < nBlock
                    cline = strip(readline(fid))
                    cline = split(cline)
                    num = length(cline)
                    for i = 1:num
                        nd = nd + 1
                        dip[nd] = parse(Float64,cline[i])
                    end
                end # while

            # slant angle
            elseif occursin("slant", cline)
                slant = zeros(nBlock)
                nd = 0
                while nd < nBlock
                    cline = strip(readline(fid))
                    cline = split(cline)
                    num = length(cline)
                    for i = 1:num
                        nd = nd + 1
                        slant[nd] = parse(Float64, cline[i])
                    end
                end # while

            end # if condType

        end # if

    end # while !eof(fid)

    close(fid)

    # convert to linear conductivity
    if condType == "isotropy"
        if modType == "log" # log10-based
            sigma = 10 .^ sigma
        end
        if resType == "resistivity"
            sigma = 1 ./ sigma
        end

    elseif condType == "anisotropy"
        if modType == "log"
            sigmax = 10 .^ sigmax
            sigmay = 10 .^ sigmay
            sigmaz = 10 .^ sigmaz
        end
        if resType == "resistivity"
            sigmax = 1 ./ sigmax
            sigmay = 1 ./ sigmay
            sigmaz = 1 ./ sigmaz
        end

    end

    ## assembly
    airDep = sum(airLayer)
    seaDep = sum(seaLayer)
    zLen = [reverse(airLayer, dims=1); seaLayer; zLen]
    nz   = length(zLen)

    if isempty(airLayer)
        airMat = zeros(0)
    else
        airMat = ones(nx*ny*nAir) * sigAir
        origin[3] = origin[3] + airDep
    end
    if isempty(seaLayer)
        seaMat = zeros(0)
    else
        seaMat = ones(nx*ny*nSea) * sigWater
    end

    if condType == "isotropy"
        sigma = vcat(airMat, seaMat, sigma)
        offsigma = zeros(0)

    elseif condType == "anisotropy"
        sigma = hcat(sigmax, sigmay, sigmaz)
        asMat = repeat(vcat(airMat, seaMat), 1, 3)
        if isempty(strike) && isempty(dip) && isempty(slant)
            offsigma = zeros(0)
            sigma = vcat(asMat, sigma)
        else
            nzb = nz - (nAir + nSea)
            aniAngle = zeros(nx*ny*nzb, 3)
            if !isempty(strike)
                aniAngle[:, 1] = strike
            end
            if !isempty(dip)
                aniAngle[:, 2] = dip
            end
            if !isempty(slant)
                aniAngle[:, 3] = slant
            end
            (sigma, offsigma) = eulerRotation!(sigma, aniAngle)
            # add airlayer and seaLayer
            sigma    = vcat(asMat, sigma)
            offsigma = vcat(zeros(size(asMat)), offsigma)
        end
    end

    # mesh statistics
    gridSize = vec([nx ny nz])
    spMat  = spzeros(0, 0)
    emMesh = EMTensorMesh(xLen, yLen, zLen, airLayer, seaLayer, gridSize, origin,
                          condType, sigma, offsigma, spMat, spMat, spMat, spMat,
                          spMat, spMat, spMat)

    return emMesh

end
