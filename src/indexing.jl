
@inline function getindex(x::SMatrix{4, 4}, i::T, j::T) where T <: Union{DNA, RNA}
    return x[trailing_zeros(i) + 1, trailing_zeros(j) + 1]
end

@inline function getindex(x::MMatrix{4, 4}, i::T, j::T) where T <: Union{DNA, RNA}
    return x[trailing_zeros(i) + 1, trailing_zeros(j) + 1]
end
