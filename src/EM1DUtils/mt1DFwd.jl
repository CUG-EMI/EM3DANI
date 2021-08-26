export mt1DAnalyticFields, mt1DAnalyticFieldsAni
export mt1DAnalyticImpAniPek, mt1DAnalyticImpAni
# export getEffSigma, rotZ

#-------------------------------------------------------------------------------
"""
`mt1DAnalyticFields` computes analytic fields for 1D layered model.
    1) time dependence: e^{iwt};
    2) wavenumber: k^2 = -iwu/p + uew^2;
    3) for Ex-Hy.
 Assuming a halfspace with the same conductivity as the bottom layer.
 Note: length(sigma) = length(zNode)-1

Input:
    freq    :: Real     - a single frequency value.
    sigma   :: Vector   - conductivity of layered model.
    zNode   :: Vector   - depth of top of each layer.
    srcType :: String   - "E" or "H", determines the top E- or H-field is source.
    fTop    :: Complex  - the given top boundary E-field or H-field value.
    compH   :: Bool     - indicates whether compute the H-field or not.

Output:
     eField :: Vector   - electric fields at the top of each layer.
     hField :: Vector   - magnetic fields at the top of each layer.

"""
function mt1DAnalyticFields(freq::T, sigma::Vector{T}, zNode::Vector{T};
                       srcType::AbstractString="E", fTop::ComplexF64=1.0+0im,
                       compH::Bool=false) where {T<:Real}

    # check if layer conductivity and depth array have the correct size
    if length(sigma) != length(zNode)-1
        error("The dimension of layer conductivity and depth mismatch!")
    end

    # physical constant
    MU0   = 4 * pi * 1e-7
    EPS0  = 8.85 * 1e-12
    omega = 2 * pi * freq
    omu   = omega * MU0

    nLayer = length(zNode)
    zLen   = diff(zNode)

    # 1st, compute the impedance at the top layer using the well-known recurrence formula for impedances.

    # wave number of all layers
    # assume halfspace at the bottom of model: [sigma; sigma[end]]
    # ka = sqrt.(MU0*EPS0*omega^2 - 1im*omu*sigma)
    ka = sqrt.(-1im*omu*[sigma; sigma[end]])


    # intrinsic impedance of the bottom layer
    ztmp = omega * MU0 / ka[end]

    for j = nLayer-1:-1:1
        k = ka[j]
	    zp = omega * MU0 / k
	    ztmp = zp * (ztmp + zp * tanh(k*zLen[j]*1im)) / (zp + ztmp * tanh(k*zLen[j]*1im))
    end
    z1 = ztmp

    # Up- and down-going electric fileds of each layer boundary
    eLayer = zeros(ComplexF64, 2, nLayer)
    hLayer = copy(eLayer)

    # Up- and down-going components of the top layer, respectively
    k1 = ka[1]

    if srcType == "E"          # E1 = c
        eLayer[1, 1] = 0.5 * fTop * (1 - omu / (z1*k1))
        eLayer[2, 1] = 0.5 * fTop * (1 + omu / (z1*k1))

    #    hLayer[1, 1] = -k1/omu * eLayer[1, 1]
    #    hLayer[2, 1] =  k1/omu * eLayer[2, 1]

    elseif srcType == "H"     # H1 = c
        hLayer[1, 1] = 0.5 * fTop * (1 - z1*k1 / omu)
        hLayer[2, 1] = 0.5 * fTop * (1 + z1*k1 / omu)

        eLayer[1, 1] = -omu/k1 * hLayer[1, 1]
        eLayer[2, 1] =  omu/k1 * hLayer[2, 1]
    end


    # 2nd, propagate the EM fields from top to bottom
    for j = 1:nLayer-1
        kr = ka[j] / ka[j+1]
        pInv = 0.5 * [1+kr 1-kr; 1-kr 1+kr]
        eUD = [exp(ka[j]*zLen[j]*1im) 0; 0 exp(-ka[j]*zLen[j]*1im)]
        eLayer[:, j+1] = pInv * eUD * eLayer[:, j]

        # check if overflow happens
        e2 = abs(eLayer[1, j+1] + eLayer[2, j+1])
        e1 = abs(eLayer[1, j]   + eLayer[2, j])
        if e2-e1>0 || isnan(e2)
            eLayer[:, j+1:end] .= 0.0
            break
        end
    end

    eField = transpose( sum(eLayer, dims=1) )

    if compH
        hLayer = [ eLayer[1:1,:] * sparse(Diagonal(-ka)) /omu;
                   eLayer[2:2,:] * sparse(Diagonal(ka))  /omu ]

        hField = transpose( sum(hLayer, dims=1) )

        return eField, hField
    end

    return eField

end




"""
`mt1DAnalyticFieldsAni` computes analytic fields for 1D general anisotropic model.
 The solution is basd on field propagation formulas (see the document by Yuguo Li).
 This function is basically wrapped from Yuguo's FORTRAN counterpart.

Input:
    freq    :: Real     - a single frequency value.
    sigma   :: Array    - nLayer*6, conductivity tensor of layered model.
    zNode   :: Vector   - nLayer+1, depth of top of each layer.
    eTop    :: Vector   - [Ex0; Ey0], the given top boundary E-field value.
    compH   :: Bool     - whether compute the H-field or not.
    timeFac :: String   - which time factor will be used: e^{iwt} ("pos",
                          by default) or e^{-iwt} ("neg")

Output:
    Ex, Ey, Ez, Hx, Hy :: Vector{Complex}
                          - EM fields for each layer (Ex, Ey, Hx, Hy at layer
                            interfaces, Ez at layer centers).

"""
function mt1DAnalyticFieldsAniNew(freq::T, sigma::Array{T}, zNode::Vector{T},
                            eTop::Vector{ComplexF64}; compH::Bool=false,
                            timeFac::AbstractString="pos") where {T<:Real}

    # check if layer conductivity and depth array have the correct size
    if size(sigma,1) != length(zNode)-1
        error("The dimension of layer conductivity and depth mismatch!")
    end

    # physical constant
    MU0   = 4 * pi * 1e-7
    omega = 2 * pi * freq
    iom   = 1im * omega * MU0    # Time dependence: e^{-iwt}

    # assume halfspace at the bottom of model
    sig    = vcat(sigma, sigma[end:end, :])
    nLayer = length(zNode)
    zLen   = diff(zNode)

    # effective azimuthal anisotropic conductivity
    sigEff = getEffSigma(sig)

    sh = (x::Complex) -> 0.5*(exp(x) - exp(-x))   # sinh()
    ch = (x::Complex) -> 0.5*(exp(x) + exp(-x))   # cosh()

    xip = zeros(ComplexF64, nLayer)
    xim = copy(xip)

    QpNeg    = zeros(Float64, nLayer)
    QmNegRec = copy(QpNeg)

    Ex = zeros(ComplexF64, nLayer)
    Ey = copy(Ex)
    Ez = zeros(ComplexF64, nLayer-1)
    Hx = copy(Ex)
    Hy = copy(Ex)

    S = zeros(ComplexF64, 4, 4, nLayer-1)

    Zxx = [];  Zxy = [];  Zyx = [];  Zyy = []

    # loop over layers to compute the impedance at the top layer and
    # the of S matrix of each layer.
    for j = nLayer:-1:1

        Axx = sigEff[j, 1]
        Ayy = sigEff[j, 2]
        Axy = sigEff[j, 3]

        ada = Axx + Ayy
        add = Axx - Ayy

        dA12 = sqrt( add^2 + 4*Axy^2 )
        add < 0 && (dA12 = -dA12)

        A1 = 0.5 * (ada + dA12)
        A2 = 0.5 * (ada - dA12)

        kp = sqrt(-iom) * sqrt(A1)
        km = sqrt(-iom) * sqrt(A2)

        # xip, xim, Qp, and Qm for each layer are saved for later use.
        xip[j] = -kp/iom
        xim[j] = -km/iom

        QpNeg[j]    = Axy == 0 ? 0 : 2*Axy/(add+dA12)         # -Qp
        QmNegRec[j] = Axy == 0 ? 0 : 0.5*(add-dA12)/Axy       # -1/Qm


        QpDerQm = QpNeg[j] * QmNegRec[j]                # Qp/Qm
        dq = 1-QpDerQm                                  # 1-Qp/Qm

        # impedance tensor at the top of the basement
        if j==nLayer

            Zxx =  (1/xip[j] - 1/xim[j]) * QmNegRec[j]/dq
            Zxy =  (1/xip[j] - QpDerQm/xim[j])/dq
            Zyx = -(1/xim[j] - QpDerQm/xip[j])/dq
            Zyy =  (1/xip[j] - 1/xim[j]) * QpNeg[j]/dq

            continue
        end

        sp = sh(kp*zLen[j])
        sm = sh(km*zLen[j])
        cp = ch(kp*zLen[j])
        cm = ch(km*zLen[j])

        Stmp = zeros(Complex{T}, 4, 4)
        Stmp[1, 1] =  cp-QpDerQm*cm
        Stmp[1, 2] = -QmNegRec[j]*(cp-cm)
        Stmp[1, 3] =  QmNegRec[j]*(sp/xip[j]-sm/xim[j])
        Stmp[1, 4] =  sp/xip[j]-QpDerQm*sm/xim[j]
        Stmp[2, 1] =  QpNeg[j]*(cp-cm)
        Stmp[2, 2] = -QpDerQm*cp+cm
        Stmp[2, 3] =  QpDerQm*sp/xip[j]-sm/xim[j]
        Stmp[2, 4] =  QpNeg[j]*(sp/xip[j]-sm/xim[j])
        Stmp[3, 1] = -QpNeg[j]*(xip[j]*sp-xim[j]*sm)
        Stmp[3, 2] =  QpDerQm*xip[j]*sp-xim[j]*sm
        Stmp[3, 3] = -QpDerQm*cp+cm
        Stmp[3, 4] = -QpNeg[j]*(cp-cm)
        Stmp[4, 1] =  xip[j]*sp-QpDerQm*xim[j]*sm
        Stmp[4, 2] = -QmNegRec[j]*(xip[j]*sp-xim[j]*sm)
        Stmp[4, 3] =  QmNegRec[j]*(cp-cm)
        Stmp[4, 4] =  cp-QpDerQm*cm

        Stmp = Stmp/dq

        S[:, :, j] = Stmp

        a11 = Stmp[3,1]*Zxx + Stmp[3,2]*Zyx + Stmp[3,3]
        a12 = Stmp[3,1]*Zxy + Stmp[3,2]*Zyy + Stmp[3,4]
        a21 = Stmp[4,1]*Zxx + Stmp[4,2]*Zyx + Stmp[4,3]
        a22 = Stmp[4,1]*Zxy + Stmp[4,2]*Zyy + Stmp[4,4]

        # adds a small value to avoid zero
        deta = a11*a22-a12*a21  + 1.0e-100

        if isinf(deta) || isnan(deta)
            b11 = 0.
            b12 = 0.
            b21 = 0.
            b22 = 0.
        else
            b11 =  a22/deta
            b12 = -a12/deta
            b21 = -a21/deta
            b22 =  a11/deta
        end

        a11 = Stmp[1,1]*Zxx + Stmp[1,2]*Zyx + Stmp[1,3]
        a12 = Stmp[1,1]*Zxy + Stmp[1,2]*Zyy + Stmp[1,4]
        a21 = Stmp[2,1]*Zxx + Stmp[2,2]*Zyx + Stmp[2,3]
        a22 = Stmp[2,1]*Zxy + Stmp[2,2]*Zyy + Stmp[2,4]

        # Z = AB
        Zxx = a11*b11 + a12*b21
        Zxy = a11*b12 + a12*b22
        Zyx = a21*b11 + a22*b21
        Zyy = a21*b12 + a22*b22

    end  # nLayer

    z = [Zxx Zxy; Zyx Zyy]

    hTop = z \ eTop

    F = vcat(eTop, hTop)

    Ex[1] = F[1]
    Ey[1] = F[2]
    Hx[1] = F[3]
    Hy[1] = F[4]


    # loop over layers again to compute the fields.
    Ft = copy(F)              # fields at the top of current layer
    Fb = zeros(size(Ft))      # fields at the bottom of current layer
    Fc = zeros(size(Ft))      # fields at the ceter of current layer

    for j=2:nLayer

        if cond(S[:, :, j-1]) > 1/eps()
            # println("Layer# $j : the S matrix tends to be ill-conditioned!")
            Fb[1:end] .= 0
        else
            Fb = S[:, :, j-1] \ Ft
        end


        Ex[j] = Fb[1]
        Ey[j] = Fb[2]
        Hx[j] = Fb[3]
        Hy[j] = Fb[4]

        # compute the S matrix for a half of current layer
        kp = -iom * xip[j-1]
        km = -iom * xim[j-1]
        QpDerQm = QpNeg[j-1]*QmNegRec[j-1]
        dq = 1-QpDerQm
        dz = 0.5 * zLen[j-1]

        sp = sh(kp*dz)
        sm = sh(km*dz)
        cp = ch(kp*dz)
        cm = ch(km*dz)

        Sh = zeros(Complex{T}, 4, 4)
        Sh[1, 1] =  cp-QpDerQm*cm
        Sh[1, 2] = -QmNegRec[j-1]*(cp-cm)
        Sh[1, 3] =  QmNegRec[j-1]*(sp/xip[j-1]-sm/xim[j-1])
        Sh[1, 4] =  sp/xip[j-1]-QpDerQm*sm/xim[j-1]
        Sh[2, 1] =  QpNeg[j-1]*(cp-cm)
        Sh[2, 2] = -QpDerQm*cp+cm
        Sh[2, 3] =  QpDerQm*sp/xip[j-1]-sm/xim[j-1]
        Sh[2, 4] =  QpNeg[j-1]*(sp/xip[j-1]-sm/xim[j-1])
        Sh[3, 1] = -QpNeg[j-1]*(xip[j-1]*sp-xim[j-1]*sm)
        Sh[3, 2] =  QpDerQm*xip[j-1]*sp-xim[j-1]*sm
        Sh[3, 3] = -QpDerQm*cp+cm
        Sh[3, 4] = -QpNeg[j-1]*(cp-cm)
        Sh[4, 1] =  xip[j-1]*sp-QpDerQm*xim[j-1]*sm
        Sh[4, 2] = -QmNegRec[j-1]*(xip[j-1]*sp-xim[j-1]*sm)
        Sh[4, 3] =  QmNegRec[j-1]*(cp-cm)
        Sh[4, 4] =  cp-QpDerQm*cm

        Sh = Sh/dq

        # there are two ways to computes fields at layer center:
        # Method 1: from bottom to center. This is deprecated because Fb
        # may be set to zero at the very beginning.
        # Fc = Sh * Fb

        # Method 2: from top to center.
        if cond(Sh) > 1/eps()
            Fc[1:end] .= 0
        else
            Fc = Sh \ Ft
        end

        # computes E-fields at layer center
        Exc = Fc[1]
        Eyc = Fc[2]
        Ez[j-1] = -(sig[j-1, 5]*Exc + sig[j-1, 6]*Eyc)/sig[j-1, 3]

        Ft = copy(Fb)

    end   # nLayer

    if timeFac == "pos"
        conj!(Ex); conj!(Ey); conj!(Ez); conj!(Hx); conj!(Hy)
    end

    if compH
        return Ex, Ey, Ez, Hx, Hy
    end

    return Ex, Ey, Ez

end  # mt1DAnalyticFieldsAni


function mt1DAnalyticFieldsAni(freq::T, sigma::Array{T}, zNode::Vector{T},
                            eTop::Vector{ComplexF64}; compH::Bool=false,
                            timeFac::AbstractString="pos") where {T<:Real}

    # check if layer conductivity and depth array have the correct size
    if size(sigma,1) != length(zNode)-1
        error("The dimension of layer conductivity and depth mismatch!")
    end

    # physical constant
    MU0   = 4 * pi * 1e-7
    omega = 2 * pi * freq
    iom   = 1im * omega * MU0    # Time dependence: e^{-iwt}

    # assume halfspace at the bottom of model
    sig    = vcat(sigma, sigma[end:end, :])
    nLayer = length(zNode)
    zLen   = diff(zNode)

    # effective azimuthal anisotropic conductivity
    sigEff = getEffSigma(sig)

    sh = (x::Complex) -> 0.5*(exp(x) - exp(-x))   # sinh()
    ch = (x::Complex) -> 0.5*(exp(x) + exp(-x))   # cosh()

    Sprod = zeros(ComplexF64, 4, 4, nLayer)
    Sprod[:, :, end] = Matrix(1.0I, 4, 4)

    M = zeros(ComplexF64, 4, 4)

    xip = zeros(ComplexF64, nLayer)
    xim = copy(xip)

    QpNeg    = zeros(Float64, nLayer)
    QmNegRec = copy(QpNeg)

    # loop over layers to get the accumulative product of S (i.e. Sprod).
    for j = nLayer:-1:1

        Axx = sigEff[j, 1]
        Ayy = sigEff[j, 2]
        Axy = sigEff[j, 3]

        ada = Axx + Ayy
        add = Axx - Ayy

        dA12 = sqrt( add^2 + 4*Axy^2 )
        add < 0 && (dA12 = -dA12)

        A1 = 0.5 * (ada + dA12)
        A2 = 0.5 * (ada - dA12)

        kp = sqrt(-iom) * sqrt(A1)
        km = sqrt(-iom) * sqrt(A2)

        # xip, xim, Qp, and Qm for each layer are saved for later use.
        xip[j] = -kp/iom
        xim[j] = -km/iom

        QpNeg[j]    = Axy == 0 ? 0 : 2*Axy/(add+dA12)         # -Qp
        QmNegRec[j] = Axy == 0 ? 0 : 0.5*(add-dA12)/Axy       # -1/Qm


        QpDerQm = QpNeg[j] * QmNegRec[j]                # Qp/Qm
        dq = 1-QpDerQm                                  # 1-Qp/Qm

        # M of the basement
        if j==nLayer

            # 1st and 3rd columns won't be used, so set them to zeros
            M = [ 0                  1     0          QmNegRec[j];
                  0           QpNeg[j]     0                    1;
                  0   -xip[j]*QpNeg[j]     0              -xim[j];
                  0             xip[j]     0   xim[j]*QmNegRec[j] ]

            continue
        end

        sp = sh(kp*zLen[j])
        sm = sh(km*zLen[j])
        cp = ch(kp*zLen[j])
        cm = ch(km*zLen[j])

        S = zeros(Complex{T}, 4, 4)
        S[1, 1] =  cp-QpDerQm*cm
        S[1, 2] = -QmNegRec[j]*(cp-cm)
        S[1, 3] =  QmNegRec[j]*(sp/xip[j]-sm/xim[j])
        S[1, 4] =  sp/xip[j]-QpDerQm*sm/xim[j]
        S[2, 1] =  QpNeg[j]*(cp-cm)
        S[2, 2] = -QpDerQm*cp+cm
        S[2, 3] =  QpDerQm*sp/xip[j]-sm/xim[j]
        S[2, 4] =  QpNeg[j]*(sp/xip[j]-sm/xim[j])
        S[3, 1] = -QpNeg[j]*(xip[j]*sp-xim[j]*sm)
        S[3, 2] =  QpDerQm*xip[j]*sp-xim[j]*sm
        S[3, 3] = -QpDerQm*cp+cm
        S[3, 4] = -QpNeg[j]*(cp-cm)
        S[4, 1] =  xip[j]*sp-QpDerQm*xim[j]*sm
        S[4, 2] = -QmNegRec[j]*(xip[j]*sp-xim[j]*sm)
        S[4, 3] =  QmNegRec[j]*(cp-cm)
        S[4, 4] =  cp-QpDerQm*cm

        S = S/dq

        Sprod[:,:,j] = S * Sprod[:,:,j+1]

    end  # nLayer


    SM = Sprod[:,:,1] * M

    # assume Ex0 & Ey0 are known, compute Cm and Dm of the bottom layer
    Ex0 = eTop[1]
    Ey0 = eTop[2]
    G = [SM[1,2]  SM[1,4]; SM[2,2]  SM[2,4]]

    #CD = G \ [Ex0; Ey0]
    CD = inv(G) * [Ex0; Ey0]
    #CD    = zeros(ComplexF64, 2, 1)
    #CD[1] = ( G[2,2]*Ex0 - G[1,2]*Ey0)/det(G)
    #CD[2] = (-G[2,1]*Ex0 + G[1,1]*Ey0)/det(G)

#     # Alternatively, one can assume Hx0 & Hy0 are known
#     Hx0 = 1.
#     Hy0 = 1.
#     G = [SM[3,2]  SM[3,4]; SM[4,2]  SM[4,4]]
#     CD = inv(G) * [Hx0; Hy0]

    # the full C of the bottom layer
    Cbot = [0; CD[1]; 0; CD[2]]

    MC = M * Cbot

    Ex = zeros(ComplexF64, nLayer)
    Ey = copy(Ex)
    Ez = zeros(ComplexF64, nLayer-1)
    Hx = copy(Ex)
    Hy = copy(Ex)

    # loop over layers again to compute the fields.
    for j=1:nLayer
        F = Sprod[:,:,j] * MC
        Ex[j] = F[1]
        Ey[j] = F[2]
        Hx[j] = F[3]
        Hy[j] = F[4]

        if j>1
            kp = -iom * xip[j-1]
            km = -iom * xim[j-1]
            QpDerQm = QpNeg[j-1]*QmNegRec[j-1]
            dq = 1-QpDerQm
            dz = 0.5 * zLen[j-1]

            sp = sh(kp*dz)
            sm = sh(km*dz)
            cp = ch(kp*dz)
            cm = ch(km*dz)

            S = zeros(Complex{T}, 4, 4)
            S[1, 1] =  cp-QpDerQm*cm
            S[1, 2] = -QmNegRec[j-1]*(cp-cm)
            S[1, 3] =  QmNegRec[j-1]*(sp/xip[j-1]-sm/xim[j-1])
            S[1, 4] =  sp/xip[j-1]-QpDerQm*sm/xim[j-1]
            S[2, 1] =  QpNeg[j-1]*(cp-cm)
            S[2, 2] = -QpDerQm*cp+cm
            S[2, 3] =  QpDerQm*sp/xip[j-1]-sm/xim[j-1]
            S[2, 4] =  QpNeg[j-1]*(sp/xip[j-1]-sm/xim[j-1])
            S[3, 1] = -QpNeg[j-1]*(xip[j-1]*sp-xim[j-1]*sm)
            S[3, 2] =  QpDerQm*xip[j-1]*sp-xim[j-1]*sm
            S[3, 3] = -QpDerQm*cp+cm
            S[3, 4] = -QpNeg[j-1]*(cp-cm)
            S[4, 1] =  xip[j-1]*sp-QpDerQm*xim[j-1]*sm
            S[4, 2] = -QmNegRec[j-1]*(xip[j-1]*sp-xim[j-1]*sm)
            S[4, 3] =  QmNegRec[j-1]*(cp-cm)
            S[4, 4] =  cp-QpDerQm*cm

            S = S/dq

            # computes E-fields at layer center
            SF = S * F
            Exc = SF[1]
            Eyc = SF[2]
            Ez[j-1] = -(sig[j-1, 5]*Exc + sig[j-1, 6]*Eyc)/sig[j-1, 3]
        end

    end   # nLayer

    if timeFac == "pos"
        conj!(Ex); conj!(Ey); conj!(Ez); conj!(Hx); conj!(Hy)
    end

    if compH
        return Ex, Ey, Ez, Hx, Hy
    end

    return Ex, Ey, Ez

end  # mt1DAnalyticFieldsAniOld


"""
`mt1DAnalyticImpAniPek` computes MT transfer function (e.g., impedance, apparent
 resistivity, phase) analytically for 1D general anisotropic model.
 The solution is basd on impedance propagation formulas (Josef Pek et al., 2002).
 This routine was basically wrapped from Pek's Fortran code 'z1adr'.

 Assumptions:
   1) The basement is a halfspace with the conductivity of the last layer.

Input:
    freqs   :: Array    - frequency array.
    sigma   :: Array    - nLayer*6, conductivity tensor of layered model.
    zNode   :: Vector   - nLayer, depth of top of each layer.
    ( Note: length(sigma) = length(zNode) )
    timeFac :: String   - which time factor will be used: e^{iwt} ("pos",
                          by default) or e^{-iwt} ("neg")

Output:
    Z, Rho, Phs :: Array - impedance, apparent resistivity, phase at the top of
                           the first layer.
"""
function mt1DAnalyticImpAniPek(freqs::Array{T}, sigma::Array{T}, zNode::Array{T};
                               timeFac::AbstractString="pos") where {T<:Real}

    # check if layer conductivity and depth array have the correct size
    if size(sigma,1) != length(zNode)
        error("The dimension of layer conductivity and depth mismatch!")
    end

    # effective azimuthal anisotropic conductivity
    sigEff = getEffSigma(sigma)

    nFreq  = length(freqs)
    nLayer = size(sigma, 1)

    h = diff(zNode)

    dfp = (x::Complex) ->  1.0 + exp(-2.0*x)
    dfm = (x::Complex) ->  1.0 - exp(-2.0*x)

    MU0 = 4*pi*1e-7


    Z   = zeros(ComplexF64, 4, nFreq)
    Rho = zeros(Float64, 4, nFreq)
    Phs = copy(Rho)

    for iFreq=1:nFreq
        omega = 2 * pi * freqs[iFreq]
        iom = 1im * omega * MU0       # Time dependence: e^{-iwt}
        k0  = sqrt(-iom)

        z = zeros(ComplexF64, 2, 2)
        Zrot = copy(z)
        bsRef = []

        for j = nLayer:-1:1

            Axx = sigEff[j, 1]
            Ayy = sigEff[j, 2]
            Axy = sigEff[j, 3]

            ada = Axx + Ayy
            add = Axx - Ayy

            dA12 = sqrt( add^2 + 4*Axy^2 )
            A1 = 0.5 * (ada + dA12)
            A2 = 0.5 * (ada - dA12)

            bs = (dA12 >= floatmin(Float64)) ? 0.5*acos(add/dA12) : 0.0

            Axy < 0 && (bs = -bs)

            kp = k0 * sqrt(A1)
            km = k0 * sqrt(A2)

            if j==nLayer
                Zrot[:, :] .= 0.0
                Zrot[1, 2] =  k0 / sqrt(A1)
                Zrot[2, 1] = -k0 / sqrt(A2)
                z = rotZ(Zrot, -bs)
                bsRef = copy(bs)
                continue
            end

            dtZbot = det(Zrot)

            if bs!=bsRef  &&  A1!=A2
                Zbot = rotZ(Zrot, bs-bsRef)
            else
                Zbot = copy(Zrot)
                bs = copy(bsRef)
            end

            dz1 = k0/sqrt(A1)
            dz2 = k0/sqrt(A2)
            ag1 = kp * h[j]
            ag2 = km * h[j]

            zdenom = dtZbot * dfm(ag1) * dfm(ag2) / (dz1*dz2)  +
                     Zbot[1,2] * dfm(ag1) * dfp(ag2) / (dz1)   -
                     Zbot[2,1] * dfp(ag1) * dfm(ag2) / (dz2)   +
                     dfp(ag1) * dfp(ag2)

            Zrot[1,1] = 4 * Zbot[1,1]*exp(-ag1-ag2)

            Zrot[1,2] = Zbot[1,2] * dfp(ag1) * dfp(ag2)            -
                        Zbot[2,1] * dfm(ag1) * dfm(ag2) * dz1/dz2  +
                        dtZbot * dfp(ag1) * dfm(ag2) /dz2          +
                        dfm(ag1) * dfp(ag2) * dz1

            Zrot[2,1] = Zbot[2,1] * dfp(ag1) * dfp(ag2)            -
                        Zbot[1,2] * dfm(ag1) * dfm(ag2) * dz2/dz1  -
                        dtZbot * dfm(ag1) * dfp(ag2) /dz1          -
                        dfp(ag1) * dfm(ag2) * dz2

            Zrot[2,2] = 4 * Zbot[2,2] * exp(-ag1-ag2)

            Zrot = Zrot / zdenom

            bsRef = copy(bs)

        end  # nLayer

        if nLayer>1
            z = bsRef==0 ? copy(Zrot) : rotZ(Zrot, -bsRef)
        end

        Zxx = z[1,1]; Zxy = z[1,2]; Zyx = z[2,1]; Zyy = z[2,2]

        rhoxx = abs(Zxx)^2 / (omega * MU0)
        rhoxy = abs(Zxy)^2 / (omega * MU0)
        rhoyx = abs(Zyx)^2 / (omega * MU0)
        rhoyy = abs(Zyy)^2 / (omega * MU0)

        phsxx = atan(imag(Zxx), real(Zxx)) * 180/pi
        phsxy = atan(imag(Zxy), real(Zxy)) * 180/pi
        phsyx = atan(imag(Zyx), real(Zyx)) * 180/pi
        phsyy = atan(imag(Zyy), real(Zyy)) * 180/pi

        Z[1:4, iFreq]   = [Zxx Zxy Zyx Zyy]
        Rho[1:4, iFreq] = [rhoxx rhoxy rhoyx rhoyy]
        Phs[1:4, iFreq] = [phsxx phsxy phsyx phsyy]

    end  # nFreq

    if timeFac == "pos"
        conj!(Z)
        Phs = -Phs
    end

    return Z, Rho, Phs

end  # mt1DAnalyticImpAniPek


"""
`mt1DAnalyticImpAni` is an alternative to `mt1DAnalyticImpAniPek`. It's also
 basd on impedance propagation formulas (see Yuguo Li's documents, 2009), which
 are slightly different from that of Josef Pek (2002).

 Assumptions, inputs and outputs are the same as that of `mt1DAnalyticImpAniPek`.
"""
function mt1DAnalyticImpAni(freqs::Array{T}, sigma::Array{T}, zNode::Array{T};
                            timeFac::AbstractString="pos") where {T<:Real}

    # check if layer conductivity and depth array have the correct size
    if size(sigma,1) != length(zNode)
        error("The dimension of layer conductivity and depth mismatch!")
    end

    # effective azimuthal anisotropic conductivity
    sigEff = getEffSigma(sigma)

    nFreq  = length(freqs)
    nLayer = size(sigma, 1)

    h = diff(zNode)

    sh = (x::Complex) -> 0.5*(exp(x) - exp(-x))   # sinh()
    ch = (x::Complex) -> 0.5*(exp(x) + exp(-x))   # cosh()

    MU0 = 4*pi*1e-7


    Z   = zeros(ComplexF64, 4, nFreq)
    Rho = zeros(Float64, 4, nFreq)
    Phs = copy(Rho)

    for iFreq=1:nFreq
        omega = 2 * pi * freqs[iFreq]
        iom = 1im * omega * MU0       # Time dependence: e^{-iwt}

        Zxx=[]; Zxy=[]; Zyx=[]; Zyy=[]
        for j = nLayer:-1:1

            Axx = sigEff[j, 1]
            Ayy = sigEff[j, 2]
            Axy = sigEff[j, 3]

            ada = Axx + Ayy
            add = Axx - Ayy

            dA12 = sqrt( add^2 + 4*Axy^2 )
            add < 0 && (dA12 = -dA12)

            A1 = 0.5 * (ada + dA12)
            A2 = 0.5 * (ada - dA12)

            kp = sqrt(-iom) * sqrt(A1)
            km = sqrt(-iom) * sqrt(A2)

            xip = -kp/iom
            xim = -km/iom

            QpNeg    = Axy == 0 ? 0 : 2*Axy/(add+dA12)         # -Qp
            QmNegRec = Axy == 0 ? 0 : 0.5*(add-dA12)/Axy       # -1/Qm

            QpDerQm = QpNeg * QmNegRec                         # Q_p/Q_m
            dq = 1-QpDerQm                                     # 1-Q_p/Q_m

            # impedance tensor at the top of the basement
            if j==nLayer
                Zxx =  (1.0/xip - 1.0/xim) * QmNegRec/dq
                Zxy =  (1.0/xip - QpDerQm/xim)/dq
                Zyx = -(1.0/xim - QpDerQm/xip)/dq
                Zyy =  (1.0/xip - 1.0/xim) * QpNeg/dq

                continue
            end

            sp = sh(kp*h[j])
            sm = sh(km*h[j])
            cp = ch(kp*h[j])
            cm = ch(km*h[j])

            S = zeros(ComplexF64, 4, 4)
            S[1, 1] =  cp-QpDerQm*cm
            S[1, 2] = -QmNegRec*(cp-cm)
            S[1, 3] =  QmNegRec*(sp/xip-sm/xim)
            S[1, 4] =  sp/xip-QpDerQm*sm/xim
            S[2, 1] =  QpNeg*(cp-cm)
            S[2, 2] = -QpDerQm*cp+cm
            S[2, 3] =  QpDerQm*sp/xip-sm/xim
            S[2, 4] =  QpNeg*(sp/xip-sm/xim)
            S[3, 1] = -QpNeg*(xip*sp-xim*sm)
            S[3, 2] =  QpDerQm*xip*sp-xim*sm
            S[3, 3] = -QpDerQm*cp+cm
            S[3, 4] = -QpNeg*(cp-cm)
            S[4, 1] =  xip*sp-QpDerQm*xim*sm
            S[4, 2] = -QmNegRec*(xip*sp-xim*sm)
            S[4, 3] =  QmNegRec*(cp-cm)
            S[4, 4] =  cp-QpDerQm*cm

            S = S/dq

            a11 = S[3,1]*Zxx + S[3,2]*Zyx + S[3,3]
            a12 = S[3,1]*Zxy + S[3,2]*Zyy + S[3,4]
            a21 = S[4,1]*Zxx + S[4,2]*Zyx + S[4,3]
            a22 = S[4,1]*Zxy + S[4,2]*Zyy + S[4,4]

            # adds a small value to avoid zero
            deta = a11*a22-a12*a21  + 1.0e-100

            b11 =  a22/deta
            b12 = -a12/deta
            b21 = -a21/deta
            b22 =  a11/deta

            a11 = S[1,1]*Zxx + S[1,2]*Zyx + S[1,3]
            a12 = S[1,1]*Zxy + S[1,2]*Zyy + S[1,4]
            a21 = S[2,1]*Zxx + S[2,2]*Zyx + S[2,3]
            a22 = S[2,1]*Zxy + S[2,2]*Zyy + S[2,4]

            # Z = AB
            Zxx = a11*b11 + a12*b21
            Zxy = a11*b12 + a12*b22
            Zyx = a21*b11 + a22*b21
            Zyy = a21*b12 + a22*b22

        end  # nLayer


        rhoxx = abs(Zxx)^2 / (omega * MU0)
        rhoxy = abs(Zxy)^2 / (omega * MU0)
        rhoyx = abs(Zyx)^2 / (omega * MU0)
        rhoyy = abs(Zyy)^2 / (omega * MU0)

        phsxx = atan(imag(Zxx), real(Zxx)) * 180/pi
        phsxy = atan(imag(Zxy), real(Zxy)) * 180/pi
        phsyx = atan(imag(Zyx), real(Zyx)) * 180/pi
        phsyy = atan(imag(Zyy), real(Zyy)) * 180/pi

        Z[1:4, iFreq]   = [Zxx Zxy Zyx Zyy]
        Rho[1:4, iFreq] = [rhoxx rhoxy rhoyx rhoyy]
        Phs[1:4, iFreq] = [phsxx phsxy phsyx phsyy]

    end  # nFreq

    if timeFac == "pos"
        conj!(Z)
        Phs = -Phs
    end

    return Z, Rho, Phs

end  # mt1DAnalyticImpAni


#
# `getEffSigma` computes effective azimuthal anisotropic conductivity
# (see Josef Pek et al., 2002).
function getEffSigma(sigma::Array{T,2}) where {T<:Real}

    sigEff = zeros(eltype(sigma), size(sigma,1), 3)
    sigEff[:,1] = sigma[:,1] - (sigma[:,5].^2) ./ sigma[:,3]
    sigEff[:,2] = sigma[:,2] - (sigma[:,6].^2) ./ sigma[:,3]
    sigEff[:,3] = sigma[:,4] - (sigma[:,5] .* sigma[:,6]) ./ sigma[:,3]

    return sigEff
end


# `rotZ` rotates the impedance ZIN by BETA (in radians) to obtain ZOUT
function rotZ(Zin::Array{Complex{T}, 2}, beta::T) where {T<:Real}

    co2 = cos(2 * beta)
    si2 = sin(2 * beta)

    sum1 = Zin[1,1] + Zin[2,2]
    sum2 = Zin[1,2] + Zin[2,1]
    dif1 = Zin[1,1] - Zin[2,2]
    dif2 = Zin[1,2] - Zin[2,1]

    Zout = copy(Zin)
    Zout[1,1] = 0.5 * ( sum1+dif1*co2 + sum2*si2)
    Zout[1,2] = 0.5 * ( dif2+sum2*co2 - dif1*si2)
    Zout[2,1] = 0.5 * (-dif2+sum2*co2 - dif1*si2)
    Zout[2,2] = 0.5 * ( sum1-dif1*co2 - sum2*si2)

    return Zout
end  # rotZ
