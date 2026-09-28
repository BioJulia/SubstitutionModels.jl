"""
    nucleotide_index(nt::NucleicAcid)

Return the position of an unambiguous nucleotide in A/C/G/T (DNA) or A/C/G/U (RNA)
order. Ambiguous nucleotides and gaps throw `ArgumentError`.
"""
function nucleotide_index(nt::NucleicAcid)
  (nt == DNA_A || nt == RNA_A) && return 1
  (nt == DNA_C || nt == RNA_C) && return 2
  (nt == DNA_G || nt == RNA_G) && return 3
  (nt == DNA_T || nt == RNA_U) && return 4
  throw(ArgumentError("Expected an unambiguous nucleotide (A, C, G, T/U), got $nt"))
end

"""
    NucleotideView(a)

Wrap a four-element vector or a 4×4 matrix with nucleotide indexing in A/C/G/T(U)
order, sharing its storage. Only one-based axes are supported. Integer indexing
and iteration follow the parent array; a single index on a matrix is linear.
Ambiguities and gaps throw `ArgumentError`. Mutation requires a mutable parent.
Use `parent` to retrieve the array or `copy` to wrap a separate copy.
"""
struct NucleotideView{T,N,A<:AbstractArray{T,N}} <: AbstractArray{T,N}
  data::A
  function NucleotideView(a::AbstractArray{T,N}) where {T,N}
    ((N == 1 && size(a) == (4,)) || (N == 2 && size(a) == (4,4))) ||
      throw(ArgumentError("Expected a four-element vector or a 4×4 matrix"))
    all(ax -> ax == Base.OneTo(4), axes(a)) ||
      throw(ArgumentError("Nucleotide views require one-based axes"))
    new{T,N,typeof(a)}(a)
  end
end

NucleotideView(a::NucleotideView) = a
Base.parent(a::NucleotideView) = a.data
Base.size(a::NucleotideView) = size(parent(a))
Base.axes(a::NucleotideView) = axes(parent(a))
Base.IndexStyle(::Type{NucleotideView{T,N,A}}) where {T,N,A} = IndexStyle(A)
Base.getindex(a::NucleotideView, inds::Vararg{Int}) = getindex(parent(a), inds...)
function Base.setindex!(a::NucleotideView, x, inds::Vararg{Int})
  setindex!(parent(a), x, inds...)
  return a
end
Base.copy(a::NucleotideView) = NucleotideView(copy(parent(a)))
Base.similar(a::NucleotideView, ::Type{T}, dims::Dims) where T = similar(parent(a), T, dims)

Base.dataids(a::NucleotideView) = Base.dataids(parent(a))
Base.unaliascopy(a::NucleotideView) = NucleotideView(Base.unaliascopy(parent(a)))
Base.to_index(a::NucleotideView, nt::NucleicAcid) = nucleotide_index(nt)
Base.checkbounds(::Type{Bool}, a::NucleotideView, nt::NucleicAcid) =
  checkbounds(Bool, parent(a), nucleotide_index(nt))
Base.checkbounds(::Type{Bool}, a::NucleotideView, i::NucleicAcid, j) =
  checkbounds(Bool, parent(a), nucleotide_index(i), j)
Base.checkbounds(::Type{Bool}, a::NucleotideView, i, j::NucleicAcid) =
  checkbounds(Bool, parent(a), i, nucleotide_index(j))
Base.checkbounds(::Type{Bool}, a::NucleotideView, i::NucleicAcid, j::NucleicAcid) =
  checkbounds(Bool, parent(a), nucleotide_index(i), nucleotide_index(j))
