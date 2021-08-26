export parsolveEM3DFwd
#-------------------------------------------------------------------------------
"""
The top-level function for solving 3D EM forward problem in parallel mode.

Input:
    emMesh   :: EMTensorMesh      - model parameters and mesh properties.
    emdata   :: EMData            - survey/data configurations.
    refModel :: RefModel          - reference model.
    lsParm   :: LinearSolverParm  - configurations of linear solver.
    pids     :: Vector            - ID of available processors.
    epField  :: Array             - primary electric fields at grid sampling points.
    fwdOnly  :: Bool              - set to true to enable the forward only mode.
    parMethod :: Int              - 1 or 2, indicates two different parallel-
                                    ization schemes.
    secondary:: Bool              - whether or not use the secondary field approach.

Output:
    if fwdOnly
        fwdResp  :: Array   - forward responses of all freqs, data components, etc.
    else
        predData :: Vector  - selected forward responses matching the data points.
        eField   :: Array   - solved electric fields (at grid sampling points).
    end
"""
function parsolveEM3DFwd(emMesh::EMTensorMesh, emData::EMData, refModel::RefModel,
             lsParm::LinearSolverParm, pids::Vector{Int}, epField::Array=zeros(0);
             fwdOnly::Bool=true, secondary::Bool=true, parMethod::Int=2)


    # Whether or not magnetic field is needed:
    #  - For MT, and CSEM in FO (forward only) mode, yes!
    #  - For CSEM not in FO mode, it depends on whether B field occurs in data type.
    # ns represents MT's two polarization modes, or CSEM's number of Tx.
    needB = false
    if typeof(emData) <: CSEMData
        probType = "csem"
        secondary = true
        ns = size(emData.txLoc, 1)
        if fwdOnly
            needB = true
        else
            for j = 1:length(emData.dataType)
                dt = emData.dataType[j]
                if occursin("B", dt)
                    needB = true
                    break
                end
            end
        end
    elseif typeof(emData) <: MTData
        probType = "mt"
        ns = 2
        needB = true
    end

    if secondary
        println("Using a secondary field approach!")
    else
        println("Using a total field approach!")
    end

    nFreq = length(emData.freqArray)
    np = length(pids)

    noEp = false
    isempty(epField) && (noEp = true)

    # parallelization Scheme 1: multiple frequencies, regardless of tranmitters.
    if typeof(lsParm) <: DirectSolverParm || ns<2 || rem(nFreq, np)==0 || parMethod==1
        println("Forward computation is parallelized over frequency!")
        if np > nFreq
            @warn("There are more processors than frequencies!")
            rmprocs(pids[nFreq+1:end])
            pids = pids[1:nFreq]
        end
        np = length(pids)

        freqIdxChn = RemoteChannel(()->Channel{Int}(nFreq))
        @async begin
            for i=1:nFreq
                put!(freqIdxChn, i)
            end
        end

        #
        if secondary && !noEp
            epRef = Array{RemoteChannel}(undef, nFreq)
            for i = 1:nFreq
                epRef[i] = initRemoteChannel(epField[:, :, i])
            end
        end

        subFwdResp = Array{Any}(undef, np)
        subeField  = Array{Any}(undef, np)
        subfreqIdx = Array{Any}(undef, np)

        @sync begin
            for i=1:np  #p in pids
                @async begin
                    if noEp
                        subfreqIdx[i], subFwdResp[i], subeField[i] =
                        remotecall_fetch(solveFwdFreq, pids[i], emMesh, emData,
                        refModel, lsParm, needB, freqIdxChn, secondary=secondary)
                    else
                        subfreqIdx[i], subFwdResp[i], subeField[i] =
                        remotecall_fetch(solveFwdFreq, pids[i], emMesh, emData,
                        refModel, lsParm, needB, freqIdxChn, epRef, secondary=secondary)
                    end
                end # @async
            end
        end # @sync

        # assemble the results
        T = eltype(subFwdResp[1])
        nRx  = size(subFwdResp[1], 1)
        nCol = size(subFwdResp[1], 2)
        fwdResp = Array{T}(undef, nRx, nCol, nFreq)

        nE = size(subeField[1], 1)
        eField = zeros(ComplexF64, nE, ns, nFreq)
        for i=1:np
            fIdx = subfreqIdx[i]
            fwdResp[:, :, fIdx] = subFwdResp[i]
            eField[:, :, fIdx]  = subeField[i]
        end


    # parallelization Scheme 2: multiple (transmitter frequency) combinations.
    else
        println("Forward computation is parallelized over frequency & transmitter!")
        nFwd  = nFreq * ns
        if np > nFwd
            @warn("There are more processors than tasks!")
            rmprocs(pids[nFwd+1:end])
            pids = pids[1:nFwd]
        end
        np = length(pids)

        freqTxIdxChn = RemoteChannel(()->Channel{Tuple}(nFwd))
        @async begin
            for i=1:nFreq
                for j=1:ns
                    put!(freqTxIdxChn, (j, i))
                end
            end
        end

        #
        if !noEp
            epRef = Array{RemoteChannel}(undef, ns, nFreq)
            for i=1:nFreq
                for j=1:ns
                    epRef[j, i] = initRemoteChannel(epField[:, j, i])
                end
            end
        end

        subeField = Array{Any}(undef, np)
        subbField = Array{Any}(undef, np)
        noEp && ( subepField = Array{Any}(undef, np) )
        subfreqTxIdx = Array{Any}(undef, np)

        @sync begin
            for i=1:np  #p in pids
                @async begin
                    if noEp
                        subfreqTxIdx[i], subeField[i], subbField[i], subepField[i] =
                        remotecall_fetch(solveFwdFreqTx, pids[i], emMesh, emData,
                        refModel, lsParm, needB, freqTxIdxChn, secondary=secondary)
                    else
                        subfreqTxIdx[i], subeField[i], subbField[i], =
                        remotecall_fetch(solveFwdFreqTx, pids[i], emMesh, emData,
                        refModel, lsParm, needB, freqTxIdxChn, epRef, secondary=secondary)
                    end
                end # @async
            end
        end # @sync

        # assemble the results
        isempty(emMesh.curlM) && getDiscreteOperators!(emMesh)

        (nF, nE) = size(emMesh.curlM)
        eField = zeros(ComplexF64, nE, ns, nFreq)
        bField = needB ? zeros(ComplexF64, nF, ns, nFreq) : zeros(ComplexF64, 0, 0, 0)
        secondary && noEp && (epField = zeros(ComplexF64, nE, ns, nFreq))

        for i=1:np
            n = length(subfreqTxIdx[i])
            for j=1:n
                sIdx, fIdx = subfreqTxIdx[i][j]
                eField[:, sIdx, fIdx]  = subeField[i][:, j]
                secondary && noEp && ( epField[:, sIdx, fIdx] = subepField[i][:, j] )
                needB && ( bField[:, sIdx, fIdx] = subbField[i][:, j] )
            end
        end

        if probType == "csem"
            # size: (nRx, nTx, nFreq)
            Ex, Ey, Ez, Bx, By, Bz = getTotalFieldsAtRxsCSEM(emMesh, emData, refModel,
                                                             eField, bField, epField)
        elseif probType == "mt"
            Ex, Ey, Ez, Bx, By, Bz = getTotalFieldsAtRxsMT(emMesh, emData, eField, bField)
            #Ex, Ey, Ez, Bx, By, Bz = getTotalFieldsAtRxsMT2(emMesh, emData, eField, bField)

        end

        if emData.phaseCon == "lag" # e^{-iwt}
            conj!(Ex); conj!(Ey); conj!(Ez); conj!(Bx); conj!(By); conj!(Bz)
            conj!(eField)
        end


        if probType == "csem"
            if needB
                fwdResp = hcat(Ex, Ey, Ez, Bx, By, Bz)
            else
                fwdResp = hcat(Ex, Ey, Ez)
            end

        elseif probType == "mt"
            fwdResp = getTransFunction(emData.freqArray, Ex, Ey, Bx, By, Bz,
                                       emData.dataType)
        end

    end  # Choosing par Scheme


    if fwdOnly
        return fwdResp, eField
    else
        predData = matchData(emData, fwdResp)
        return predData, eField
    end

end


#
# solve forward problem looping over frequencies.
#
function solveFwdFreq(emMesh::EMTensorMesh, emData::EMData, refModel::RefModel,
                      lsParm::LinearSolverParm, needB::Bool,
                      freqIdxChn::RemoteChannel, epRef::Array=zeros(0);
                      secondary::Bool=true)

    if typeof(emData) <: CSEMData
        probType = "csem"
        ns = size(emData.txLoc, 1)
    elseif typeof(emData) <: MTData
        probType = "mt"
        ns = 2
    end

    noEp = false
    isempty(epRef) && (noEp = true)

    if isempty(emMesh.curlM)
        @printf("%-40s", "Forming discretized operators:")
        t = @elapsed emMesh = getDiscreteOperators!(emMesh)
        @printf("%8.2f %s\n", t, "seconds.")
    end


    # extract things
    curlM = emMesh.curlM
    gradM = copy(emMesh.gradM)    # gradM may be changed later
    aveFC = emMesh.aveFC
    aveEC = emMesh.aveEC
    aveNC = emMesh.aveNC
    volM  = emMesh.volM
    sigma = emMesh.sigma
    gridSize = emMesh.gridSize
    condType = emMesh.condType
    phaseCon = emData.phaseCon

    nx, ny, nz = gridSize
    nNode = (nx+1)*(ny+1)*(nz+1)

    MU0 = 4 * pi * 1e-7
    nGrid = prod(gridSize)
    muMat = ones(nGrid) * MU0

    # setup mass matrices
    MmuF = aveFC' * (volM * (1 ./ muMat))
    MmuF = sparse(Diagonal(MmuF))

    # get discretized reference/background model
    secondary && ( sigBG = get3dBGModel(emMesh, refModel) )

    # get mass matrices for normal and abnormal/secondary conductivities
    if condType == "isotropy"
        println("This is an isotropic case! ")
        Msig  = aveEC' * (volM * sigma)
        Msig  = sparse(Diagonal(Msig))

        if secondary
            MsigS = aveEC' * (volM * (sigma - sigBG))
            MsigS = sparse(Diagonal(MsigS))
        end

    elseif condType == "anisotropy"
        println("This is an anisotropic case! ")
        Msig  = getAniMassMatrix(emMesh, sigma)
        secondary && ( MsigS = getAniMassMatrix(emMesh, sigma - sigBG) )

    end

    # prepare coefficient matrix
    Ke_re = (curlM' * MmuF) * curlM
    Ke_im = Msig

    # probType == "mt" && ( ii, io = getBoundIndex(gridSize) )
    if probType == "mt"
        ii, io = getBoundIndex(gridSize)
    end

    if typeof(lsParm) <: DirectSolverParm
        println("Linear Solver: $(uppercase( String(lsParm.solverName) ))")

    elseif typeof(lsParm) <: IterativeSolverParm
        println("Linear Solver: $(uppercase( String(lsParm.iterMethod) ))")

        # iteration log file
        cid = @sprintf("%03d", myid())
        iterlogfile = "id" * cid * "_iter.log"
        iterLogID = open(iterlogfile, "w")

        # prepare stuff for preconditioner
        if lsParm.prec == :aphi
            MmuN = aveNC' * (volM * (1 ./ muMat))
            MmuN = sparse(Diagonal(MmuN))
            STBa = gradM * MmuN * gradM'    # -∇μ⁻¹∇̇̇.

            if probType == "mt"
                # get the inner node gradient.
                idx3D  = reshape(collect(1:nNode), nx+1, ny+1, nz+1)
                iiND   = idx3D[2:end-1, 2:end-1, 2:end-1][:]
                gradM = gradM[ii, iiND]
                STBa  = STBa[ii, ii]
                Msig  = Msig[ii, ii]
            end
        end

    end


    # predefined electric and magnetic fields to be solved, and other things
    # to be output
    (nF, nE) = size(curlM)
    eField  = zeros(ComplexF64, 0)
    bField  = zeros(ComplexF64, 0)
    secondary && ( epField = zeros(ComplexF64, 0) )
    freqIdxArray = Array{Int}(undef, 0)

    nFreq = 0
    # loop over frequency
    while true

        !isready(freqIdxChn) && break

        eSingle  = zeros(ComplexF64, nE, ns)

        freqIdx = take!(freqIdxChn)
        append!(freqIdxArray, freqIdx)
        nFreq += 1

        freqArray = [ emData.freqArray[freqIdx] ]
        omega = 2 * pi * freqArray

        if secondary
            epSingle = zeros(ComplexF64, nE, ns)
            if noEp
                @printf("%-40s", "Calculating primary electric fields:")
                t = @elapsed epSingle[:, :] = compPrimaryElectricField(emMesh, emData, refModel, freqIdx)
                @printf("%8.2f %s\n", t, "seconds.")
            else
                epSingle = take!(epRef[freqIdx])
            end
        end


        # For MT, boundary fields are non-zero
        if probType == "mt"
            @printf("%-40s", "Calculating boundary fields:")
            # compute boundary total fields for two polarization modes
            t = @elapsed begin
                bcXY = getBoundaryMT(freqArray, emMesh, 1)
                bcYX = getBoundaryMT(freqArray, emMesh, 2)
            end
            @printf("%8.2f %s\n", t, "seconds.")

            eSingle[io, 1] = bcXY[:, 1]
            eSingle[io, 2] = bcYX[:, 1]

            if secondary
                # get boundary primary fields
                bcXYp = similar(bcXY)
                bcYXp = similar(bcYX)
                bcXYp[:, 1] = epSingle[io, 1]
                bcYXp[:, 1] = epSingle[io, 2]

                # get boundary secondary fields
                bcXY -= bcXYp
                bcYX -= bcYXp
            end

        end

        iome = 1im * omega[1]

        Ke = Ke_re + iome * Ke_im

        secondary && ( rhs = -iome * MsigS * epSingle )

        if probType == "mt"
            Keb = Ke[ii, io]
            Ke  = Ke[ii, ii]

            if secondary
                rhs = rhs[ii, :] - Keb * hcat(bcXY[:, 1], bcYX[:, 1])
            else
                rhs = - Keb * hcat(bcXY[:, 1], bcYX[:, 1])
            end
        end


        println("Solving the forward linear system of the $freqIdx th frequency ...")

        # using cholesky decomposition to slove the sytem of linear equation
        if typeof(lsParm) <: DirectSolverParm

            t = @elapsed etmp, = directSolver(Ke, rhs, lsParm)
            @printf("%s %8.2f %s\n", "   -", t, " seconds in total.")

        # using Krylov method to slove the sytem of linear equation
        elseif typeof(lsParm) <: IterativeSolverParm

            if lsParm.prec == :aphi        # SSOR preconditioner based on A-phi system
                Aap = [Ke + STBa              iome * Msig * gradM;
                       iome * gradM' * Msig   iome * gradM' * Msig * gradM]

                Map(x)  = tril(Aap) \ ( diag(Aap) .* (triu(Aap) \ x) )
                #@time P1 = [spunit(size(Ke,2));  gradM']
                #P1 = transpose( [spunit(size(Ke,2))  gradM] )
                P2 = [spunit(size(Ke,1))   gradM]
                P1 = transpose(P2)
                PC(x)  = P2*(Map(P1*x))

            elseif lsParm.prec == :ssor   # standard SSOR preconditioner
                #D  = diag(Ke)
                #PC(x) = tril(Ke)\(diag(Ke) .* (triu(Ke)\x))
                # PC(x) = D.\x  # Jacobi

            end

            @printf(iterLogID,"%s %3d\n","Freq No.:", freqIdx)
            t = @elapsed etmp = iterativeSolver(Ke, rhs, lsParm, iterLogID, PC)
            @printf("%s %8.2f %s\n", "   -", t, " seconds in total.")

        end


        # For MT, we perform the Rx field interpolation in terms of total fields
        if probType == "mt"
            if secondary
                eSingle[ii, :] = etmp + epSingle[ii, :]
            else
                eSingle[ii, :] = etmp
            end

        else
            eSingle[:, :]  = etmp
        end


        if needB
            bSingle = -1 / (omega[1] * 1im) * (curlM * eSingle)
        end

        if nFreq == 1
            eField  = copy(eSingle)
            secondary && ( epField = copy(epSingle) )
            needB && ( bField = copy(bSingle) )
        else
            eField  = hcat(eField, eSingle)
            secondary && ( epField = hcat(epField, epSingle) )
            needB && ( bField = hcat(bField, bSingle) )
        end


    end  # while true

    typeof(lsParm) <: IterativeSolverParm  &&  close(iterLogID)


    eField  = reshape(eField, nE, ns, nFreq)
    secondary && ( epField = reshape(epField, nE, ns, nFreq) )
    needB && ( bField = reshape(bField, nF, ns, nFreq) )


    if probType == "csem"
        # size: (nRx, nTx, nFreq)
        Ex, Ey, Ez, Bx, By, Bz = getTotalFieldsAtRxsCSEM(emMesh, emData, refModel,
                                                         eField, bField, epField, freqIdxArray)
    elseif probType == "mt"
        Ex, Ey, Ez, Bx, By, Bz = getTotalFieldsAtRxsMT(emMesh, emData, eField, bField, freqIdxArray)
        #Ex, Ey, Ez, Bx, By, Bz = getTotalFieldsAtRxsMT2(emMesh, emData, eField, bField)
    end

    if phaseCon == "lag" # e^{-iwt}
        conj!(Ex); conj!(Ey); conj!(Ez); conj!(Bx); conj!(By); conj!(Bz)
        conj!(eField)
    end


    if probType == "csem"
        if needB
            fwdResp = hcat(Ex, Ey, Ez, Bx, By, Bz)
        else
            fwdResp = hcat(Ex, Ey, Ez)
        end

    elseif probType == "mt"
        freqs = emData.freqArray[freqIdxArray]
        fwdResp = getTransFunction(freqs, Ex, Ey, Bx, By, Bz, emData.dataType)
    end

    return freqIdxArray, fwdResp, eField

end


#
# solve forward problem looping over (transmitter, frequency) combinations.
#
function solveFwdFreqTx(emMesh::EMTensorMesh, emData::EMData, refModel::RefModel,
                        lsParm::IterativeSolverParm, needB::Bool,
                        freqTxIdxChn::RemoteChannel, epRef::Array=zeros(0);
                        secondary::Bool=true)

    if typeof(emData) <: CSEMData
        probType = "csem"
    elseif typeof(emData) <: MTData
        probType = "mt"
    end

    noEp = false
    isempty(epRef) && (noEp = true)

    if isempty(emMesh.curlM)
        @printf("%-40s", "Forming discretized operators:")
        t = @elapsed emMesh = getDiscreteOperators!(emMesh)
        @printf("%8.2f %s\n", t, "seconds.")
    end


    # extract things
    curlM = emMesh.curlM
    gradM = copy(emMesh.gradM)    # gradM may be changed later
    aveFC = emMesh.aveFC
    aveEC = emMesh.aveEC
    aveNC = emMesh.aveNC
    volM  = emMesh.volM
    sigma = emMesh.sigma
    gridSize = emMesh.gridSize
    condType = emMesh.condType

    nx, ny, nz = gridSize
    nNode = (nx+1)*(ny+1)*(nz+1)

    MU0 = 4 * pi * 1e-7
    nGrid = prod(gridSize)
    muMat = ones(nGrid) * MU0

    # setup mass matrices
    MmuF = aveFC' * (volM * (1 ./ muMat))
    MmuF = sparse(Diagonal(MmuF))

    # get discretized reference/background model
    secondary && ( sigBG = get3dBGModel(emMesh, refModel) )

    # get mass matrices for normal and abnormal/secondary conductivities
    if condType == "isotropy"
        println("This is an isotropic case! ")
        Msig  = aveEC' * (volM * sigma)
        Msig  = sparse(Diagonal(Msig))

        if secondary
            MsigS = aveEC' * (volM * (sigma - sigBG))
            MsigS = sparse(Diagonal(MsigS))
        end

    elseif condType == "anisotropy"
        println("This is an anisotropic case! ")
        Msig  = getAniMassMatrix(emMesh, sigma)
        secondary && ( MsigS = getAniMassMatrix(emMesh, sigma - sigBG) )

    end

    # prepare coefficient matrix
    Ke_re = (curlM' * MmuF) * curlM
    Ke_im = Msig

    # probType == "mt" && ( ii, io = getBoundIndex(gridSize) )
    if probType == "mt"
        ii, io = getBoundIndex(gridSize)
    end

    # use iterative solver
    println("Linear Solver: $(uppercase( String(lsParm.iterMethod) ))")

    # iteration log file
    cid = @sprintf("%03d", myid())
    iterlogfile = "id" * cid * "_iter.log"
    iterLogID = open(iterlogfile, "w")

    # prepare stuff for preconditioner
    if lsParm.prec == :aphi
        MmuN = aveNC' * (volM * (1 ./ muMat))
        MmuN = sparse(Diagonal(MmuN))
        STBa = gradM * MmuN * gradM'    # -∇μ⁻¹∇̇̇.

        if probType == "mt"
            # get the inner node gradient.
            idx3D  = reshape(collect(1:nNode), nx+1, ny+1, nz+1)
            iiND   = idx3D[2:end-1, 2:end-1, 2:end-1][:]
            gradM = gradM[ii, iiND]
            STBa  = STBa[ii, ii]
            Msig  = Msig[ii, ii]
        end
    end


    # predefined electric and magnetic fields to be solved, and other things
    # to be output
    (nF, nE) = size(curlM)
    eField  = zeros(ComplexF64, 0)
    bField  = zeros(ComplexF64, 0)
    epField = zeros(ComplexF64, 0)
    #freqTxIdx = Array{Tuple}(undef, 0)
    freqTxIdx = []

    nFwd = 0
    # loop over (txIdx, freqIdx) pairs
    while true

        !isready(freqTxIdxChn) && break

        eSingle  = zeros(ComplexF64, nE)
        bSingle  = zeros(ComplexF64, nF)

        ftIdx = take!(freqTxIdxChn)
        txIdx, freqIdx = ftIdx
        append!(freqTxIdx, [ftIdx])

        nFwd += 1
        #sleep(rand(1:5))

        if secondary
            epSingle = zeros(ComplexF64, nE)
            if noEp
                @printf("%-40s", "Calculating primary electric fields:")
                t = @elapsed epSingle[:] = compPrimaryEfieldFreqTx(emMesh, emData, refModel, [ftIdx])
                @printf("%8.2f %s\n", t, "seconds.")
            else
                epSingle[:] = take!(epRef[txIdx, freqIdx])
            end
        end


        # For MT, boundary fields are non-zero
        if probType == "mt"
            @printf("%-40s", "Calculating boundary fields:")

            # bc = zeros(ComplexF64, length(io), nFwd)
            t = @elapsed bc = getBoundaryMT(emData.freqArray[[freqIdx]], emMesh, txIdx)
            @printf("%8.2f %s\n", t, "seconds.")

            eSingle[io] = bc[:]

            if secondary
                # get boundary primary fields
                bcp = similar(bc)
                bcp[:] = epSingle[io]

                # get boundary secondary fields
                bc -= bcp
            end

        end


        freq = emData.freqArray[freqIdx]
        omega = 2 * pi * freq

        iome = 1im * omega

        Ke = Ke_re + iome * Ke_im

        secondary && ( rhs = -iome * MsigS * epSingle )

        if probType == "mt"
            Keb = Ke[ii, io]
            Ke  = Ke[ii, ii]

            if secondary
                rhs = rhs[ii] - Keb * bc
            else
                rhs = - Keb * bc
            end
        end

        println("Solving the forward linear system at freqID=$freqIdx, txID=$txIdx.")
        if lsParm.prec == :aphi        # SSOR preconditioner based on A-phi system
            Aap = [Ke + STBa              iome * Msig * gradM;
                   iome * gradM' * Msig   iome * gradM' * Msig * gradM]

            Map(x)  = tril(Aap) \ ( diag(Aap) .* (triu(Aap) \ x) )
            #@time P1 = [spunit(size(Ke,2));  gradM']
            #P1 = transpose( [spunit(size(Ke,2))  gradM] )
            P2 = [spunit(size(Ke,1))   gradM]
            P1 = transpose(P2)
            PC(x)  = P2*(Map(P1*x))

        elseif lsParm.prec == :ssor   # standard SSOR preconditioner
            #D  = diag(Ke)
            #PC(x) = tril(Ke)\(diag(Ke) .* (triu(Ke)\x))
            # PC(x) = D.\x  # Jacobi

        end

        @printf(iterLogID,"%s %3d\n","Freq No.:", freqIdx)
        t = @elapsed etmp = iterativeSolver(Ke, rhs, lsParm, iterLogID, PC)
        @printf("%s %8.2f %s\n", "   -", t, " seconds in total.")

        # For MT, we perform the Rx field interpolation in terms of total fields
        if probType == "mt"
            if secondary
                eSingle[ii] = etmp + epSingle[ii]
            else
                eSingle[ii] = etmp
            end

        else
            eSingle = etmp
        end


        append!(eField, eSingle)

        if needB
            bSingle = -1 / iome * (curlM * eSingle)
            append!(bField, bSingle)
        end

        secondary && noEp && append!(epField, epSingle)

    end  # while true

    close(iterLogID)

    eField = reshape(eField, nE, nFwd)
    needB && ( bField  = reshape(bField,  nF, nFwd) )
    secondary && noEp && ( epField = reshape(epField, nE, nFwd) )

    return freqTxIdx, eField, bField, epField

end
