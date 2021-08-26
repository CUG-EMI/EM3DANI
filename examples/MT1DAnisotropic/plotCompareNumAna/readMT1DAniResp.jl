export MT1DResp, readMT1DAniResp

mutable struct MT1DResp{T<:Real}
    freqs   :: Vector{T}
    Imp     :: Array{T}
    appRho  :: Array{T}
end

#-------------------------------------------------------------------------------
"""
Read MT 1D anisotropic forward responses from file.

Input:
    respfile :: String   - the name of the response file.

Output:
    resp     :: MT1DResp - forward response.

"""
function readMT1DAniResp(respfile::String)

    if isfile(respfile)
        fid = open(respfile, "r")
    else
        error("$(respfile) does not exist, please try again.")
    end

    freqs = [];    Imp = [];   appRho = []

    nFreq = 0;

    while !eof(fid)

        cline = strip(readline(fid))

        # ignore all comments: empty line, or line preceded with #
        while isempty(cline) || cline[1] == '#'
            cline = strip(readline(fid))
        end

        # frequencies
        if occursin("Frequencies", cline)
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

        # Impedance
        elseif occursin("Impedance", cline)

            nd  = 0
            Imp = zeros(Float64, nFreq, 8)

            while nd < nFreq
                cline = strip(readline(fid))
                while isempty(cline) || cline[1] == '#'
                    cline = strip(readline(fid))
                end
                cline = split(cline)
                nd = nd + 1

                for j=1:8
                    Imp[nd, j] = parse(Float64, cline[j])
                end

            end

        # Apparent resistivity
        elseif occursin("Apparent resistivity", cline)

            nd  = 0
            appRho = zeros(Float64, nFreq, 8)

            while nd < nFreq
                cline = strip(readline(fid))
                while isempty(cline) || cline[1] == '#'
                    cline = strip(readline(fid))
                end
                cline = split(cline)
                nd = nd + 1

                for j=1:8
                    appRho[nd, j] = parse(Float64, cline[j])
                end

            end

        end

    end

    close(fid)


    resp = MT1DResp(freqs, Imp, appRho)

    return resp

end
