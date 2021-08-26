# Construct mass matrix for anisotropic conductivity model.
function getAniMassMatrix(emMesh::EMTensorMesh, sigma::Array{Float64, 2})

    # mesh info
    gridSize = emMesh.gridSize
    aveEC    = emMesh.aveEC

    # anisotropic conductivity
    offsigma = emMesh.offsigma
    volM     = emMesh.volM

    # number of edges
    nGrid = prod(gridSize)
    nx, ny, nz = gridSize
    nEx = nx * (ny+1) * (nz+1)
    nEy = (nx+1) * ny * (nz+1)
    nEz = (nx+1) * (ny+1) * nz
    #
    aveEx = aveEC[:, 1:nEx]
    aveEy = aveEC[:, nEx+1:nEx+nEy]
    aveEz = aveEC[:, nEx+nEy+1:end]
    aveEx = transpose(aveEx) * volM
    aveEy = transpose(aveEy) * volM
    aveEz = transpose(aveEz) * volM
    Mxxsig = aveEx * sigma[:, 1]
    Myysig = aveEy * sigma[:, 2]
    Mzzsig = aveEz * sigma[:, 3]
    Msig  = sparse(Diagonal(vcat(Mxxsig, Myysig, Mzzsig)))
    if isempty(offsigma)
        printstyled("tri-axis anisotropy! \n", color=:cyan)
    else
        printstyled("general anisotropy! \n", color=:cyan)
        sigxy = offsigma[:, 1]
        sigxz = offsigma[:, 2]
        sigyz = offsigma[:, 3]
        MsigOD = getOffDiagMassMatrix(gridSize, volM, sigxy, sigxz, sigyz)
        Msig += MsigOD
    end

    return Msig

end



# Get the part arising from off-diagonal entries of the conductiviy tensor.
function getOffDiagMassMatrix(gridSize::Vector{Int}, Vol::SparseMatrixCSC,
            sigxy::Vector{T}, sigxz::Vector{T}, sigyz::Vector{T}) where {T<:Real}


    nx, ny, nz = gridSize
    nEx = nx*(ny+1)*(nz+1)
    nEy = (nx+1)*ny*(nz+1)
    nEz = (nx+1)*(ny+1)*nz

    # define some utility anonymous functions
    # av   = (n::Int) -> spdiagm((0.5*ones(n), 0.5*ones(n)), (0, 1), n, n+1)
    oeye = (n::Int) -> vcat(spzeros(1,n), spunit(n))     # size: (n+1,n)
    eyeo = (n::Int) -> vcat(spunit(n), spzeros(1,n))

    # Ex -> Ey
    a_xy1 = kron(spunit(nz+1), kron(oeye(ny), av(nx)))  # size: (nEx,nEy)
    a_xy2 = kron(spunit(nz+1), kron(eyeo(ny), av(nx)))
    A_xy1 = [spzeros(nEx,nEx)    a_xy1    spzeros(nEx,nEz)]
    A_xy2 = [spzeros(nEx,nEx)    a_xy2    spzeros(nEx,nEz)]

    #
    avnx2t = copy(transpose(2 * av(nx)))
    avny2t = copy(transpose(2 * av(ny)))
    avnz2t = copy(transpose(2 * av(nz)))

    c2E_xy1 = kron(avnz2t, kron(oeye(ny), spunit(nx)))  # size: (nEx,nCell)
    c2E_xy2 = kron(avnz2t, kron(eyeo(ny), spunit(nx)))
    Msigxy1 = sdiag(c2E_xy1 * 0.25 * Vol * sigxy)
    Msigxy2 = sdiag(c2E_xy2 * 0.25 * Vol * sigxy)

    Axy = Msigxy1 * A_xy1 + Msigxy2 * A_xy2  # size: (nEx,nE)


    # Ex -> Ez
    a_xz1 = kron(oeye(nz), kron(spunit(ny+1), av(nx)))  # size: (nEx,nEz)
    a_xz2 = kron(eyeo(nz), kron(spunit(ny+1), av(nx)))
    A_xz1 = [spzeros(nEx,nEx)    spzeros(nEx,nEy)    a_xz1]
    A_xz2 = [spzeros(nEx,nEx)    spzeros(nEx,nEy)    a_xz2]

    c2E_xz1 = kron(oeye(nz), kron(avny2t, spunit(nx))) # size: (nEx,nCel)
    c2E_xz2 = kron(eyeo(nz), kron(avny2t, spunit(nx)))
    Msigxz1 = sdiag(c2E_xz1 * 0.25 * Vol * sigxz)
    Msigxz2 = sdiag(c2E_xz2 * 0.25 * Vol * sigxz)

    Axz = Msigxz1 * A_xz1 + Msigxz2 * A_xz2  # size: (nEx,nE)

    AX = Axy + Axz


    # Ey -> Ex
    a_yx1 = kron(spunit(nz+1), kron(av(ny), oeye(nx)))  # size: (nEy,nEx)
    a_yx2 = kron(spunit(nz+1), kron(av(ny), eyeo(nx)))
    A_yx1 = [a_yx1    spzeros(nEy,nEy)    spzeros(nEy,nEz)]
    A_yx2 = [a_yx2    spzeros(nEy,nEy)    spzeros(nEy,nEz)]

    c2E_yx1 = kron(avnz2t, kron(spunit(ny), oeye(nx)))  # size: (nEy,nCell)
    c2E_yx2 = kron(avnz2t, kron(spunit(ny), eyeo(nx)))
    Msigyx1 = sdiag(c2E_yx1 * 0.25 * Vol * sigxy)
    Msigyx2 = sdiag(c2E_yx2 * 0.25 * Vol * sigxy)

    Ayx = Msigyx1 * A_yx1 + Msigyx2 * A_yx2  # size: (nEy,nE)


    # Ey -> Ez
    a_yz1 = kron(oeye(nz), kron(av(ny), spunit(nx+1)))  # size: (nEy,nEz)
    a_yz2 = kron(eyeo(nz), kron(av(ny), spunit(nx+1)))
    A_yz1 = [spzeros(nEy,nEx)    spzeros(nEy,nEy)    a_yz1]
    A_yz2 = [spzeros(nEy,nEx)    spzeros(nEy,nEy)    a_yz2]

    c2E_yz1 = kron(oeye(nz), kron(spunit(ny), avnx2t)) # size: (nEy,nCell)
    c2E_yz2 = kron(eyeo(nz), kron(spunit(ny), avnx2t))
    Msigyz1 = sdiag(c2E_yz1 * 0.25 * Vol * sigyz)
    Msigyz2 = sdiag(c2E_yz2 * 0.25 * Vol * sigyz)

    Ayz = Msigyz1 * A_yz1 + Msigyz2 * A_yz2 # size: (nEx,nE)

    AY = Ayx + Ayz


    # Ez -> Ex
    a_zx1 = kron(av(nz), kron(spunit(ny+1), oeye(nx)))  # size: (nEz,nEx)
    a_zx2 = kron(av(nz), kron(spunit(ny+1), eyeo(nx)))
    A_zx1 = [a_zx1    spzeros(nEz,nEy)    spzeros(nEz,nEz)]
    A_zx2 = [a_zx2    spzeros(nEz,nEy)    spzeros(nEz,nEz)]

    c2E_zx1 = kron(spunit(nz), kron(avny2t, oeye(nx)))  # size: (nEz,nCell)
    c2E_zx2 = kron(spunit(nz), kron(avny2t, eyeo(nx)))
    Msigzx1 = sdiag(c2E_zx1 * 0.25 * Vol * sigxz)
    Msigzx2 = sdiag(c2E_zx2 * 0.25 * Vol * sigxz)

    Azx = Msigzx1 * A_zx1 + Msigzx2 * A_zx2  # size: (nEz,nE)


    # Ez -> Ey
    a_zy1 = kron(av(nz), kron(oeye(ny), spunit(nx+1)))  # size: (nEz,nEy)
    a_zy2 = kron(av(nz), kron(eyeo(ny), spunit(nx+1)))
    A_zy1 = [spzeros(nEz,nEx)    a_zy1    spzeros(nEz,nEz)]
    A_zy2 = [spzeros(nEz,nEx)    a_zy2    spzeros(nEz,nEz)]

    c2E_zy1 = kron(spunit(nz), kron(oeye(ny), avnx2t))  # size: (nEz,nCell)
    c2E_zy2 = kron(spunit(nz), kron(eyeo(ny), avnx2t))
    Msigzy1 = sdiag(c2E_zy1 * 0.25 * Vol * sigyz)
    Msigzy2 = sdiag(c2E_zy2 * 0.25 * Vol * sigyz)

    Azy = Msigzy1 * A_zy1 + Msigzy2 * A_zy2  # size: (nEz,nE)

    AZ = Azx + Azy

    return vcat(AX, AY, AZ)

end
