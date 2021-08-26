export readEdgeFields
#-------------------------------------------------------------------------------

"""
Read grid edge vector variables (e.g. E-field) from file.

Input:
    fieldfile :: String   - the name of the file for input.

Output:
    eField    :: Array    - grid edge vector variables.
    gridSize  :: Vector   - cell numbers in x,y,z direction.

"""
function readEdgeFields(fieldfile::String)

    if isfile(fieldfile)
        fid = open(fieldfile, "r")
    else
        error("$fieldfile does not exist, please try again.")
    end

	# Read some basic information
	cline = strip(readline(fid))
	tmp   = split(cline)
	nFreq = parse(Int, tmp[end])

	cline = strip(readline(fid))
	tmp   = split(cline)
	nTx   = parse(Int, tmp[end])

	cline = strip(readline(fid))
	tmp   = split(cline)
	gridSize = zeros(Int, 3)
	gridSize[1] = parse(Int, tmp[end-2])
	gridSize[2] = parse(Int, tmp[end-1])
	gridSize[3] = parse(Int, tmp[end])

	nx, ny, nz = gridSize
	nEx = nx * (ny+1) * (nz+1)
	nEy = (nx+1) * ny * (nz+1)
	nEz = (nx+1) * (ny+1) * nz
	nE  = nEx+nEy+nEz

    eField = zeros(ComplexF64, nE, nTx, nFreq)

  # Loop over frequencies and transmitters to read the fields.
	cline = readline(fid)
    for i=1:nFreq
		  cline = readline(fid)

          for j=1:nTx
			    cline = readline(fid)

			    cline = readline(fid)
			    nd = 0
			    tmpR = zeros(nEx)
			    while nd < nEx
			        cline = strip(readline(fid));  cline = split(cline)
                    for jj = 1:length(cline)
                        nd = nd + 1
                        tmpR[nd] = parse(Float64, cline[jj])
                    end
                end

			    cline = readline(fid)
			    nd = 0
			    tmpI = zeros(nEx)
			    while nd < nEx
			        cline = strip(readline(fid));  cline = split(cline)
                    for jj = 1:length(cline)
                        nd = nd + 1
                        tmpI[nd] = parse(Float64, cline[jj])
                    end
                end

			    eField[1:nEx,j,i] = tmpR + tmpI * 1im


			    cline = readline(fid)
			    nd = 0
			    tmpR = zeros(nEy)
			    while nd < nEy
			        cline = strip(readline(fid));  cline = split(cline)
                    for jj = 1:length(cline)
                        nd = nd + 1
                        tmpR[nd] = parse(Float64, cline[jj])
                    end
                end

			    cline = readline(fid)
			    nd = 0
			    tmpI = zeros(nEy)
			    while nd < nEy
			        cline = strip(readline(fid));  cline = split(cline)
                    for jj = 1:length(cline)
                        nd = nd + 1
                        tmpI[nd] = parse(Float64, cline[jj])
                    end
                end

			    eField[nEx+1:nEx+nEy,j,i] = tmpR + tmpI * 1im


			    cline = readline(fid)
			    nd = 0
			    tmpR = zeros(nEz)
			    while nd < nEz
			        cline = strip(readline(fid));  cline = split(cline)
                    for jj = 1:length(cline)
                        nd = nd + 1
                        tmpR[nd] = parse(Float64, cline[jj])
                    end
                end

			    cline = readline(fid)
			    nd = 0
			    tmpI = zeros(nEz)
			    while nd < nEz
			        cline = strip(readline(fid));  cline = split(cline)
                    for jj = 1:length(cline)
                        nd = nd + 1
                        tmpI[nd] = parse(Float64, cline[jj])
                    end
                end

			    eField[nEx+nEy+1:nEx+nEy+nEz,j,i] = tmpR + tmpI * 1im


        end  # over nTx
    end   # over nFreq

	close(fid)

	return eField, gridSize

end
