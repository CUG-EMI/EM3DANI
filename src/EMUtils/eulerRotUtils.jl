export eulerRotation!

# Perform Euler's rotation to transform the principal conductivity to a tensor.
function eulerRotation!(sigma::Array{T}, aniAngle::Array{T}) where {T<:Real}

    #
    if size(sigma) != size(aniAngle)
        error("Dimension is not consistent between conductivity and angles.")
    end
    # anisotropic conductivity
    nGrid  = size(sigma, 1)
    strike = aniAngle[:, 1]
    dip    = aniAngle[:, 2]
    slant  = aniAngle[:, 3]

    # rotations
    sinS = sind.(strike)
    cosS = cosd.(strike)
    sinD = sind.(dip)
    cosD = cosd.(dip)
    sinL = sind.(slant)
    cosL = cosd.(slant)

    #
    sigxx = zeros(Float64, nGrid)
    sigyy = zeros(Float64, nGrid)
    sigzz = zeros(Float64, nGrid)
    sigxy = zeros(Float64, nGrid)
    sigxz = zeros(Float64, nGrid)
    sigyz = zeros(Float64, nGrid)

    @inbounds for j=1:nGrid
        rzS = [cosS[j]  -sinS[j]  0;
               sinS[j]   cosS[j]  0;
                 0         0      1]
        rxD = [1     0          0;
               0   cosD[j]  -sinD[j];
               0   sinD[j]   cosD[j]]
        rzL = [cosL[j]  -sinL[j]  0;
               sinL[j]   cosL[j]  0;
                 0         0      1]
        R = rzS * rxD * rzL
        sigAni = spdiagm(0 => sigma[j, :])
        sigAni = R * sigAni * R'
        #
        sigAni[abs.(sigAni) .< eps()] .= 0.0
        sigxx[j] = sigAni[1, 1]
        sigxy[j] = sigAni[2, 1]
        sigxz[j] = sigAni[3, 1]
        sigyy[j] = sigAni[2, 2]
        sigyz[j] = sigAni[3, 2]
        sigzz[j] = sigAni[3, 3]
    end

    sigma = hcat(sigxx, sigyy, sigzz)
    offsigma = hcat(sigxy, sigxz, sigyz)

    return sigma, offsigma

end
