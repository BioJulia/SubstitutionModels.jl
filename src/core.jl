abstract type SubstitutionModel end


const SM = SubstitutionModel


"""
`NucleicAcidSubstitutionModel` is an abstract type that contains all models
describing a substitution process impacting biological sequences of `DNA` or
`RNA` with continous time Markov models.
"""
abstract type NucleicAcidSubstitutionModel <: SubstitutionModel end


const NASM = NucleicAcidSubstitutionModel


const Qmatrix = SMatrix{4, 4, Float64}


const Pmatrix = SMatrix{4, 4, Float64}


function Base.convert(::Type{T}, θ::F...; safe::Bool=true) where {T <: NASM, F <: Float64}
  return T(θ..., safe=safe)
end


function Base.convert(::Type{T}, θ_vec::A; safe::Bool=true) where {T <: NASM, A <: AbstractArray}
  return T(θ_vec, safe=safe)
end


function Base.convert(::Type{T}, θ_vec::AbstractArray, π_vec::AbstractArray; safe::Bool=true) where T <: NASM
  return T(θ_vec, π_vec, safe=safe)
end


function Base.show(io::IO, mod::NASM)
  print(io, nameof(typeof(mod)), '(')
  for i in 1:fieldcount(typeof(mod))
    i > 1 && print(io, ", ")
    show(io, getfield(mod, i))
  end
  print(io, ')')
  return nothing
end
