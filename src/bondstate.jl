struct ExtMatrix{T<:Number}#<:AbstractMatrix{T} where {T<:Number}
    M::Matrix{T}
end

function Base.getindex(A::ExtMatrix{T},i,j)::T where {T<:Number}
    (H,L) = size(A.M)
    if i<=0||j<=0||i>=H+1||j>=L+1
        return one(T)
    else
        return A.M[i,j]
    end 
end
1+1
function Base.setindex!(A::ExtMatrix{T},v,i,j) where {T<:Number}
    if 1<=i<=size(A.M,1) && 1<=j<=size(A.M,2)
        A.M[i,j] = v
    else
        @argcheck v ≈ 1 "cannot assign value $v at inded $((i,j)) to matrix of size $(size(A.M))"
    end
end
Base.size(A) = size(A.M)
Base.length(A) = length(A.M)

mutable struct Bondstate{T<:Number}
    H::Int
    L::Int
    h::ExtMatrix{T}
    v::ExtMatrix{T}
    d::ExtMatrix{T}
    R::T
    
    function Bondstate(H::Int,L::Int,h::Matrix{T},v::Matrix{T},d::Matrix{T}) where {T}
        @argcheck size(h) == (H,L-1)
        @argcheck size(v) == (H-1,L)
        @argcheck size(d) == (H-1,L-1)
        return new{T}(H,L,ExtMatrix{T}(h),ExtMatrix{T}(v),ExtMatrix{T}(d),1) 
    end
end


