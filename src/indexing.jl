const NucleicAcidlookup8 = [0x01, 0x02, 0x00, 0x03, 0x00, 0x00, 0x00, 0x04]

function getindex{T <: Union{DNA, RNA}}(x::SMatrix{4, 4}, i::T, j::T)
  return x[NucleicAcidlookup8[reinterpret(UInt8, i)], NucleicAcidlookup8[reinterpret(UInt8, j)]]
end

function getindex{T <: Union{DNA, RNA}}(x::MMatrix{4, 4}, i::T, j::T)
  return x[NucleicAcidlookup8[reinterpret(UInt8, i)], NucleicAcidlookup8[reinterpret(UInt8, j)]]
end

function setindex!{T <: Union{DNA, RNA}}(x::MMatrix{4,4}, v, i::T, j::T)
  return setindex!(x, v, NucleicAcidlookup8[reinterpret(UInt8, i)], NucleicAcidlookup8[reinterpret(UInt8, j)])
end
