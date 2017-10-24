
@inline function Base.checkbounds(x::SMatrix{4, 4}, i::T, j::T) where T <: NucleicAcid
    if iscertain(i) && iscertain(j)
        return true
    end
    throw(BoundsError(x, i, j))
end

@inline function getindex(x::SMatrix{4, 4}, i::T, j::T) where T <: NucleicAcid
    @boundscheck checkbounds(x, i, j)
    @inbounds return x[trailing_zeros(i) + 1, trailing_zeros(j) + 1]
end

@inline function getindex(x::MMatrix{4, 4}, i::T, j::T) where T <: NucleicAcid
    @boundscheck checkbounds(x, i, j)
    @inbounds return x[trailing_zeros(i) + 1, trailing_zeros(j) + 1]
end
