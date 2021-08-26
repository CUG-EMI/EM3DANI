export solveEM3DFwd
#-------------------------------------------------------------------------------
"""
The top-level function for solving 3D EM forward problem.

Input:
    emMesh   :: EMTensorMesh      - model parameters and mesh properties.
    emdata   :: EMData            - survey/data configurations.
    refModel :: RefModel          - reference model.
    epField  :: Array             - primary electric fields at grid sampling points.
    lsParm   :: LinearSolverParm  - configurations of linear solver.
    fwdOnly  :: Bool              - set to true to enable the forward only mode.
    secondary:: Bool              - whether or not use the secondary field approach.

Output:
    if fwdOnly
        fwdResp  :: Array   - forward responses of all freqs, data components, etc.
    else
        predData :: Vector  - selected forward responses matching the data points.
        eField   :: Array   - solved electric fields (at grid sampling points).
    end
"""
function solveEM3DFwd(emMesh::EMTensorMesh, emData::EMData, refModel::RefModel,
                      lsParm::LinearSolverParm, epField::Array=zeros(0);
                      fwdOnly::Bool=true, secondary::Bool=true)


    if typeof(emData) <: CSEMData
        probType = "csem"
        secondary = true
    elseif typeof(emData) <: MTData
        probType = "mt"
    end

    if secondary
        println("Using a secondary field approach!")
    else
        println("Using a total field approach!")
    end

    if secondary && isempty(epField)
        probType == "mt"    && @printf("%-40s", "Calculating primary electric fields:")
        probType == "csem"  && println("Calculating primary electric fields ...")
        t = @elapsed epField = compPrimaryElectricField(emMesh, emData, refModel)
        probType == "mt"    && @printf("%8.2f %s\n", t, "seconds.")
    end

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

    # frequencies
    freqArray = emData.freqArray
    omega = 2 * pi * freqArray
    nFreq = length(omega)


    # Whether magnetic field is needed or not:
    #  - For MT, and CSEM in FO (forward only) mode, yes!
    #  - For CSEM not in FO mode, it depends on whether B field occurs in data type.
    # ns represents MT's two polarization modes, or CSEM's number of Tx.
    needB = false
    if probType == "mt"
        ns = 2
        needB = true

    elseif probType == "csem"
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
    end


    # predefined electric and magnetic fields to be solved
    (nF, nE) = size(curlM)
    eField = zeros(ComplexF64, nE, ns, nFreq)
    bField = needB ? zeros(ComplexF64, nF, ns, nFreq) : zeros(ComplexF64, 0, 0, 0)


    # prepare coefficient matrix
    Ke_re = (curlM' * MmuF) * curlM
    Ke_im = Msig


    # For MT, boundary fields are non-zero
    if probType == "mt"
        @printf("%-40s", "Calculating boundary fields:")
        # compute boundary total fields for two polarization modes
        t = @elapsed begin
            ii, io = getBoundIndex(gridSize)
            bcXY = getBoundaryMT(freqArray, emMesh, 1)
            bcYX = getBoundaryMT(freqArray, emMesh, 2)
        end
        @printf("%8.2f %s\n", t, "seconds.")

        eField[io, 1, :] = bcXY[:, :]
        eField[io, 2, :] = bcYX[:, :]

        if secondary
            # get boundary primary fields
            bcXYp = similar(bcXY)
            bcYXp = similar(bcYX)
            bcXYp[:, :] = epField[io, 1, :]
            bcYXp[:, :] = epField[io, 2, :]

            # get boundary secondary fields
            bcXY -= bcXYp
            bcYX -= bcYXp
        end

    end



    if typeof(lsParm) <: DirectSolverParm
        println("Linear Solver: $(uppercase( String(lsParm.solverName) ))")

        # predefined forward decomposization factors for direct solver
        if !lsParm.saveFac
            fwdDecomp = Array{Any}(undef, nFreq)
        elseif lsParm.solverName == :mumps
            fwdDecomp = Array{MUMPSfactorization}(undef, nFreq)
        elseif lsParm.solverName == :mklpardiso
            fwdDecomp = Array{MKLPardisoSolver}(undef, nFreq)
        end

    elseif typeof(lsParm) <: IterativeSolverParm
        println("Linear Solver: $(uppercase( String(lsParm.iterMethod) ))")

        # iteration log file
        iterLogID = open("iter.log", "w")

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


    # loop over frequency
    for i = 1:nFreq
        iome = 1im * omega[i]

        Ke = Ke_re + iome * Ke_im

        secondary && ( rhs = -iome * MsigS * epField[:, :, i] )

        if probType == "mt"
            Keb = Ke[ii, io]
            Ke  = Ke[ii, ii]

            if secondary
                rhs = rhs[ii, :] - Keb * hcat(bcXY[:, i], bcYX[:, i])
            else
                rhs = - Keb * hcat(bcXY[:, i], bcYX[:, i])
            end

        end


        println("Solving the forward linear system of the $i / $nFreq th frequency ...")

        # using cholesky decomposition to slove the sytem of linear equation
        if typeof(lsParm) <: DirectSolverParm

            t = @elapsed etmp, fwdDecomp[i] = directSolver(Ke, rhs, lsParm)
            @printf("%s %8.2f %s\n", "   -", t, " seconds in total.")
            # destroyMUMPS(Ainv)

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

            @printf(iterLogID,"%s %3d\n","Freq No.:", i)
            t = @elapsed etmp = iterativeSolver(Ke, rhs, lsParm, iterLogID, PC)
            @printf("%s %8.2f %s\n", "   -", t, " seconds in total.")

        end


        # For MT, we perform the Rx field interpolation in terms of total fields
        if probType == "mt"
            if secondary
                eField[ii, :, i] = etmp + epField[ii, :, i]
            else
                eField[ii, :, i] = etmp
            end

        else
            eField[:, :, i]  = etmp
        end


        if needB
            tmp = view(eField, 1:nE, 1:ns, i)
            bField[:, :, i] = -1 / (omega[i] * 1im) * (curlM * tmp)
        end

    end  # nFreq

    typeof(lsParm) <: IterativeSolverParm  &&  close(iterLogID)


    if probType == "csem"
        # size: (nRx, nTx, nFreq)
        Ex, Ey, Ez, Bx, By, Bz = getTotalFieldsAtRxsCSEM(emMesh, emData, refModel,
                                                         eField, bField, epField)
    elseif probType == "mt"
        Ex, Ey, Ez, Bx, By, Bz = getTotalFieldsAtRxsMT(emMesh, emData, eField, bField)
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
        fwdResp = getTransFunction(freqArray, Ex, Ey, Bx, By, Bz, emData.dataType)
    end


    if fwdOnly
        return fwdResp, eField
    else
        predData = matchData(emData, fwdResp)
        return predData, eField
    end

end
