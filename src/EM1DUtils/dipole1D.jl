#------------------------------------------------------------------------------
#
# dipole1D defines interface to frotran code Dipole1D (developed by Kerry key),
# which allows us to call functions to compute quasi-analytic solution of 1D
# CSEM problem.
#
# (C) Peng Ronghua, 16 Aug, 2015
#
#------------------------------------------------------------------------------
export comp1DCSEM

const dipole1dLibPath = abspath(joinpath(splitdir(Base.source_path())[1], ".",
                        "deps", "lib", "Dipole1D"))

"""
    `comp1DCSEM(txLoc, rxLoc, freqArray, sigma, depth1D)`

Calculate quasi-analytic electric fields of 1D CSEM problem due to electric dipole.

# Inputs:
 - `txLoc`     =::Array[nsx5], transmitter location
 - `rxLoc`     =::Array[nrx3], receiver location
 - `freqArray` =::Array, transmission frequencies
 - `sigma`     =::Vector, conductivities of layered model
 - `depth1D`   =::Vector, depths of top of layered model

# Outputs:
 - `ex1d`  =::Array, x-component of electric fields
 - `ey1d`  =::Array, y-component of electric fields
 - `ez1d`  =::Array, z-component of electric fields

"""
function comp1DCSEM(txLoc::Array{T}, rxLoc::Array{T}, freqArray::Array{T},
                    sigma::Vector{T}, depth1D::Vector{T}) where{T<:Real}

    # transmitter and receivers
    nTx = size(txLoc, 1)
    nRx = size(rxLoc, 1)
    nFreq  = length(freqArray)
    nLayer = length(sigma)

    if length(sigma) != length(depth1D)
        error("comp1DCSEM: Conductivity and depth must be the same length!")
    end

    # predicted field component
    ex1d = zeros(ComplexF64, nRx, nTx, nFreq)
    ey1d = zeros(ComplexF64, nRx, nTx, nFreq)
    ez1d = zeros(ComplexF64, nRx, nTx, nFreq)

    #
    ccall((:calldipole1d_, dipole1dLibPath), Nothing, (Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                    Ptr{Float64}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64},
                    Ptr{Int64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}),
                    txLoc, rxLoc, freqArray, sigma, depth1D, Ref(nTx), Ref(nRx), Ref(nFreq),
                    Ref(nLayer), ex1d, ey1d, ez1d)

    return ex1d, ey1d, ez1d

end # comp1DCSEM

"""
    `comp1DCSEM(txLoc, rxLoc, freqArray, sigma, depth1D, compFlag)`

Calculate specfic component of EM field due to electric dipole.

# Inputs:
 - `txLoc`     =::Array[nsx5], transmitter location
 - `rxLoc`     =::Array[nrx3], receiver location
 - `freqArray` =::Array, transmission frequencies
 - `sigma`     =::Vector, conductivities of layered model
 - `depth1D`   =::Vector, depths of top of layered model
 - `compFlag`  =::Integer, component identifier

# Outputs:
 - `e1d` =::Array, computed field

"""
function comp1DCSEM(txLoc::Array{T,2}, rxLoc::Array{T,2},
                    freqArray::Array{T}, sigma::Vector{T},
                    depth1D::Vector{T}, compFlag::Integer) where{T<:Real}

    # transmitter and receivers
    nTx = size(txLoc, 1)
    nRx = size(rxLoc, 1)
    nFreq  = length(freqArray)
    nLayer = length(sigma)

    if length(sigma) != length(depth1D)
        error("comp1DCSEM: Conductivity and depth must be the same length!")
    end

    # predicted field component
    e1d = zeros(ComplexF64, nRx, nTx, nFreq)

    #
    ccall((:calldipole1deb_, dipole1dLibPath), Nothing, (Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                    Ptr{Float64}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64},
                    Ptr{Int64}, Ptr{Int64}, Ptr{ComplexF64}), txLoc, rxLoc, freqArray,
                    sigma, depth1D, Ref(nTx), Ref(nRx), Ref(nFreq), Ref(nLayer), Ref(compFlag), e1d)

    return e1d

end # comp1DCSEM


"""
    `comp1DCSEM(txLoc, rxLoc, dipLen, nIntPts, freqArray, sigma, depth1D, compFlag)`

Calculate specific component of EM field due to finite direct wire source.

# Inputs:
 - `txLoc`   =::Array[nsx5], transmitter location
 - `rxLoc`   =::Array[nrx3], receiver location
 - `dipLen`  =::Float64, length of finite dipole source
 - `nIntPts` =::Integer, Number of points to use for Gauss quadrature
                integration for finite dipole
 - `freqArray` =::Array, transmission frequencies
 - `sigma`     =::Array, conductivities of layered model
 - `depth1D`   =::Array, depth of top of layered model
 - `compFlag`  =::Integer, component identifier

# Outputs:
 - `e1d` =::Array, computed EM field

"""
function comp1DCSEM(txLoc::Array{T}, rxLoc::Array{T}, dipLen::T,
                    nIntPts::Integer, freqArray::Array{T},
                    sigma::Vector{T}, depth1D::Vector{T},
                    compFlag::Integer) where{T<:Real}

    # transmitter and receivers
    nTx = size(txLoc, 1)
    nRx = size(rxLoc, 1)
    nFreq  = length(freqArray)
    nLayer = length(sigma)

    if length(sigma) != length(depth1D)
        error("comp1DCSEM: Conductivity and depth must be the same length!")
    end

    # predicted field component
    e1d = zeros(ComplexF64, nRx, nTx, nFreq)

    #
    ccall((:calldipole1dfinite_, dipole1dLibPath), Nothing, (Ptr{Float64}, Ptr{Float64},
                        Ptr{Float64}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64},
                        Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64},
                        Ptr{Int64}, Ptr{Int64}, Ptr{ComplexF64}), txLoc, rxLoc,
                        Ref(dipLen), Ref(nIntPts), freqArray, sigma, depth1D, Ref(nTx), Ref(nRx),
                        Ref(nFreq), Ref(nLayer), Ref(compFlag), e1d)

    return e1d

end # comp1DCSEM
