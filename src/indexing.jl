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

function indexing_depwarn(method::Symbol)
  Base.depwarn("Indexing arbitrary AbstractArrays with DNA/RNA symbols is deprecated " *
               "and will be removed in SubstitutionModels 0.6.0. " *
               "Use array[nucleotide_index(nt)] (convert both indices for matrices).", method)
end

function Base.checkbounds(a::AbstractArray, i::NucleicAcid)
  indexing_depwarn(:checkbounds)
  checkbounds(a, nucleotide_index(i))
end

function Base.checkbounds(a::AbstractArray, i::T, j::T) where T <: NucleicAcid
  indexing_depwarn(:checkbounds)
  checkbounds(a, nucleotide_index(i), nucleotide_index(j))
end

function Base.getindex(a::AbstractArray, i::NucleicAcid)
  indexing_depwarn(:getindex)
  return a[nucleotide_index(i)]
end

function Base.getindex(a::AbstractArray, i::T, j::T) where T <: NucleicAcid
  indexing_depwarn(:getindex)
  return a[nucleotide_index(i), nucleotide_index(j)]
end

function Base.setindex!(a::AbstractArray, x, i::NucleicAcid)
  indexing_depwarn(:setindex!)
  return setindex!(a, x, nucleotide_index(i))
end

function Base.setindex!(a::AbstractArray, x, i::T, j::T) where T <: NucleicAcid
  indexing_depwarn(:setindex!)
  return setindex!(a, x, nucleotide_index(i), nucleotide_index(j))
end
