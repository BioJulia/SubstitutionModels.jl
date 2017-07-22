const NucleicAcidlookup = [0x01, 0x02, 0x00, 0x03, 0x00, 0x00, 0x00, 0x04]

function getindex(x::Array{Any, 1}, i::NucleicAcid)
  return x[NucleicAcidlookup[reinterpret(UInt8, i)]]
end

function getindex{T <: Union{DNA, RNA}}(x::Array{Any, 2}, i::T, j::T)
  return x[NucleicAcidlookup[reinterpret(UInt8, i)], NucleicAcidlookup[reinterpret(UInt8, j)]]
end
