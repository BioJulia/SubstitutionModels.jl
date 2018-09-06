@inline function checkbounds(a::AbstractArray, i::T, j::T) where T <: NucleicAcid
  if iscertain(i) & iscertain(j)
    return true
  end
  checkbounds(a, trailing_zeros(i) + 1, trailing_zeros(j) + 1)
end


@inline function checkbounds(a::AbstractArray, i::NucleicAcid)
  if iscertain(i)
      return true
  end
  checkbounds(a, trailing_zeros(i) + 1)
end


@inline function getindex(a::AbstractArray, i::T, j::T) where T <: NucleicAcid
  @boundscheck checkbounds(a, i, j)
  @inbounds return a[trailing_zeros(i) + 1, trailing_zeros(j) + 1]
end


@inline function getindex(a::AbstractArray, i::NucleicAcid)
  @boundscheck checkbounds(a, i)
  @inbounds return a[trailing_zeros(i) + 1]
end


@inline function setindex!(a::AbstractArray, x, i::T, j::T) where T <: NucleicAcid
  @boundscheck checkbounds(a, i, j)
  @inbounds return setindex!(a, x, trailing_zeros(i) + 1, trailing_zeros(j) + 1)
end


@inline function setindex!(a::AbstractArray, x, i::NucleicAcid)
  @boundscheck checkbounds(a, i)
  @inbounds return setindex!(a, x, trailing_zeros(i) + 1)
end
