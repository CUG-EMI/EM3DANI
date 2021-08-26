
function readEM1DAniParam(paramfile::String)

    if isfile(paramfile)
        fid = open(paramfile, "r")
    else
        error("$(paramfile) does not exist, please try again.")
    end

    # Initialization
    h = []
    sigma = []
    aniAngle = []
    freqs = []

    while !eof(fid)

        cline = strip(readline(fid))

        # ignore empty lines and comments (preceded with #)
        isempty(cline)  && continue
        cline[1] == '#' && continue

        # block at x-axis
        if occursin("Layers", cline)
            tmp = split(cline)
            nLayer  = parse(Int, tmp[end])

            readline(fid)  # skip one line
            nd  = 0
            h = zeros(nLayer)
            sigma    = zeros(nLayer, 3)
            aniAngle = zeros(nLayer, 3)

            while nd < nLayer
                cline = strip(readline(fid))
                cline = split(cline)
                nd = nd + 1
                h[nd] = parse(Float64, cline[1])
                sigma[nd, 1] = parse(Float64, cline[2])
                sigma[nd, 2] = parse(Float64, cline[3])
                sigma[nd, 3] = parse(Float64, cline[4])
                aniAngle[nd, 1] = parse(Float64, cline[5])
                aniAngle[nd, 2] = parse(Float64, cline[6])
                aniAngle[nd, 3] = parse(Float64, cline[7])
            end

        # block at y-axis
        elseif occursin("Frequencies", cline)
            tmp = split(cline)
            nFreq  = parse(Int, tmp[end])
            nd  = 0
            freqs = zeros(nFreq)
            while nd < nFreq
                cline = strip(readline(fid))
                while isempty(cline) || cline[1] == '#'
                    cline = strip(readline(fid))
                end
                #cline = split(cline)
                nd = nd + 1
                freqs[nd] = parse(Float64, cline)
            end

        end # if

    end # while !eof(fid)

    close(fid)


    # from thickness to depth
    h = [0; h[1:end-1]]
    zNode = cumsum(h) * 1000

    # # from resistivity to conductivity tensor
    sigma    = 1 ./ sigma
    sigma, offsigma = eulerRotation!(sigma, aniAngle)
    sigma = hcat(sigma, offsigma)

    return freqs, sigma, zNode

end
