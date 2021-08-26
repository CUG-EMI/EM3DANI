#------------------------------------------------------------------------------
# Define some routines for linear interpolation.
#------------------------------------------------------------------------------
export linearInterp, locateNearest, locateSegment1D
export linearInterpMat, bilinearInterpMat, trilinearInterpMat


# Get 1D linear interpolation points and weights for a single point.
function linearInterp(point::T, x::Vector{T}) where {T<:Real}

    val, ind = findmin(abs.(point .- x))
    # location of point
    if point - x[ind] > 0
    # point on the right
        indL = ind
        indR = ind + 1
    else
    # point on the left
        indL = ind - 1
        indR = ind
    end

    # ensure interpolation points within the bound
    n = length(x)
    indL = maximum([minimum([indL, n]), 1])  # 1<=indL<=n
    indR = maximum([minimum([indR, n]), 1])  # 1<=indR<=n

    # interpolation weights
    if indL != indR
        xLen = x[indR] - x[indL]
        wL = 1 - (point - x[indL]) / xLen
        wR = 1 - (x[indR] - point) / xLen
    else
        wL = 0.5
        wR = 0.5
    end

    return indL, indR, wL, wR

end


# Get 1D linear interpolation points and weights and use them to construct the
# interpolation matrix.
function linearInterpMat(point::Vector{T}, x::Vector{T}, getMat=true) where {T<:Real}

    npts  = length(point)
    nNode = length(x)

    inds    = Array{Any}(undef, npts)
    weights = zeros(Float64, npts, 2)
    interpMat = spzeros(Float64, nNode, npts)

    for i = 1:npts
        indL, indR, wL, wR = linearInterp(point[i], x)
        inds[i] = (indL,indR)
        weights[i,:] = [wL wR]
        interpMat[:,i] = sparsevec([indL;indR], [wL;wR], nNode)
    end

    !getMat && ( return inds, weights )

    return interpMat

end


# Get bilinear interpolation points and weights and use them to construct the
# interpolation matrix.
function bilinearInterpMat(point::Array{T,2}, x::Vector{T}, y::Vector{T},
                           getMat=true) where {T<:Real}

    npts = size(point, 1)
    nx = length(x);  ny = length(y)
    nNode = nx*ny

    inds    = Array{Any}(undef, npts)
    weights = zeros(Float64, npts, 4)
    interpMat = spzeros(Float64, nNode, npts)

    for i = 1:npts
        indxL, indxR, wxL, wxR = linearInterp(point[i, 1], x)
        indyL, indyR, wyL, wyR = linearInterp(point[i, 2], y)
        inds[i] = [(indxL, indyL), (indxR, indyL),
                   (indxL, indyR), (indxR, indyR)]

        w1 = wxL * wyL
        w2 = wxR * wyL
        w3 = wxL * wyR
        w4 = wxR * wyR

        w = [w1 w2 w3 w4]
        weights[i,:] = w

        idx = zeros(Int, 4)
        for j = 1:4
            idx[j] = nx * (inds[i][j][2]-1) + inds[i][j][1]
        end
        interpMat[:,i] = sparsevec(idx, vec(w), nNode)

    end

    !getMat && ( return inds, weights )

    return interpMat

end


# Get trilinear interpolation points and weights and use them to construct the
# interpolation matrix.
function trilinearInterpMat(point::Array{T, 2}, x::Vector{T}, y::Vector{T},
                            z::Vector{T}, getMat=true) where {T<:Real}

    npts = size(point,1)
    nx = length(x);  ny = length(y);  nz = length(z)
    nNode = nx*ny*nz

    inds    = Array{Any}(undef, npts)
    weights = zeros(Float64, npts, 8)
    interpMat = spzeros(Float64, nNode, npts)

    for i = 1:npts
        indxL, indxR, wxL, wxR = linearInterp(point[i, 1], x)
        indyL, indyR, wyL, wyR = linearInterp(point[i, 2], y)
        indzL, indzR, wzL, wzR = linearInterp(point[i, 3], z)
        inds[i] = [(indxL, indyL, indzL),
                   (indxR, indyL, indzL),
                   (indxL, indyR, indzL),
                   (indxR, indyR, indzL),
                   (indxL, indyL, indzR),
                   (indxR, indyL, indzR),
                   (indxL, indyR, indzR),
                   (indxR, indyR, indzR)]

        w1 = wxL * wyL * wzL
        w2 = wxR * wyL * wzL
        w3 = wxL * wyR * wzL
        w4 = wxR * wyR * wzL
        w5 = wxL * wyL * wzR
        w6 = wxR * wyL * wzR
        w7 = wxL * wyR * wzR
        w8 = wxR * wyR * wzR

        w = [w1 w2 w3 w4 w5 w6 w7 w8]
        weights[i,:] = w

        idx = zeros(Int, 8)
        for j = 1:8
            idx[j] = nx * ny * (inds[i][j][3]-1) + nx * (inds[i][j][2]-1) + inds[i][j][1]
        end
        interpMat[:,i] = sparsevec(idx, vec(w), nNode)

    end

    !getMat && ( return inds, weights )

    return interpMat

end


# Find the location of points to the nearest coordinate location (node, cell-
# center, etc.).
function locateNearest(point::Array{T,2}, xLoc::Vector{T}, yLoc::Vector{T},
                       zLoc::Vector{T}) where {T<:Real}

    location = ones(Int,3)

    # x location
    val,ind = findmin(abs.(point[1] .- xLoc))
    location[1] = ind

    # y location
    val,ind = findmin(abs.(point[2] .- yLoc))
    location[2] = ind

    # z location
    val,ind = findmin(abs.(point[3] .- zLoc))
    location[3] = ind

    return location

end


# Determine the segment a point reside in. Note that if the point is at the
# boundary, assume it's in the left segment.
function locateSegment1D(xp::T, x::Vector{T}) where {T<:Real}

    n = length(x)
    val, ind = findmin(abs.(x .- xp))
    if xp <= x[ind]; ind -= 1; end
    ind = maximum([minimum([ind, n-1]), 1])

end
