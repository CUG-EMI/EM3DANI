#------------------------------------------------------------------------------
# Define necessary routines for numerical discretization.
#------------------------------------------------------------------------------
import Base: kron

export kron3
export ddx, ddxC2N, av, avcn
export spdiag, spdiags
export sdiag, sdInv, ddxCellGradBC
export ndgrid, meshgrid
export spunit


# Construct sparse unit matrix, replace the built-in function speye in old versions.
function spunit(n::Integer)
    return sparse(1.0I, n, n)
end


#
function spdiag((x1,x2), (d1,d2), m, n)
    I, J, V = SparseArrays.spdiagm_internal(d1 => x1, d2 => x2)
    return sparse(I, J, V, m, n)
end


# Form sparse diagonal matrix
function spdiags(v::AbstractArray)
    return sparse(Diagonal(v))
end


# Form sparse diagonal matrix, an alternative way
function sdiag(v::AbstractVector)
    return spdiagm(0 => v)
end


# Compute inverse matrix of a sparse diagonal matrix
function sdInv(sdMat::SparseMatrixCSC)
    M = diag(sdMat)
    return sdiag(1 ./ M)
end


# Form 1D difference sparse matrix from node to center
function ddx(n::Integer)
    return spdiag((-ones(n), ones(n)), (0,1), n, n+1)
end


# Form 1D difference sparse matrix from grid center to node
function ddxC2N(n::Integer)
    dC2N = spdiag((-ones(n), ones(n)), (-1,0), n+1, n)
    dC2N[end, end] = 1.0
    return dC2N
end


# Form 1D averaging matrix from node to cell-center
function av(n::Integer)
    return spdiag((0.5*ones(n), 0.5*ones(n)), (0,1), n, n+1)
end


# Form 1D averaging matrix from cell-center to node
function avcn(n::Integer)
    avn = spdiag((0.5*ones(n), 0.5*ones(n)), (-1,0), n+1, n)
    avn[1,1]     = 1.0
    avn[end,end] = 1.0
    return avn
end


# Compute cell-center 1D gradient operator with boundary condition
function ddxCellGradBC(n::Integer, BC::String)

    GradBC = spdiag((-ones(n), ones(n)), (-1,0), n+1, n)
    BC     = lowercase(BC)

    # first side
    if BC == "neumann"
        GradBC[1,1] = 0.
        # second side
        GradBC[end,end]  = 0.
    elseif BC == "dirichlet"
        GradBC[1,1] = -2.
        # second side
        GradBC[end,end] = 2.
    end

    return GradBC

end



# Compute kronecker tensor products for three matrices
function kron3(A::T, B::T, C::T) where {T<:SparseMatrixCSC}
    return kron(A, kron(B, C))
end


#
ndgrid(v::AbstractVector) = copy(v)

function ndgrid_fill(a, v, s, snext)
	for j = 1:length(a)
		a[j] = v[div(rem(j-1, snext), s)+1]
	end
end


function ndgrid(vs::AbstractVector{T}...) where{T}
	n = length(vs)
	sz = map(length, vs)
	out = ntuple(i->Array{T}(undef, sz), n)
	s = 1
	for i=1:n
		a = out[i]::Array
		v = vs[i]
		snext = s*size(a,i)
		ndgrid_fill(a, v, s, snext)
		s = snext
	end
	out
end


#
meshgrid(v::AbstractVector) = meshgrid(v, v)
function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T}) where{T}
	m, n = length(vy), length(vx)
	vx = reshape(vx, 1, n)
	vy = reshape(vy, m, 1)
	(repeat(vx, m, 1), repeat(vy, 1, n))
end


function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T},
                     vz::AbstractVector{T}) where{T}

	m, n, o = length(vy), length(vx), length(vz)
	vx = reshape(vx, 1, n, 1)
	vy = reshape(vy, m, 1, 1)
	vz = reshape(vz, 1, 1, o)
	om = ones(Int, m)
	on = ones(Int, n)
	oo = ones(Int, o)
	(vx[om, :, oo], vy[:, on, oo], vz[om, on, :])
end
