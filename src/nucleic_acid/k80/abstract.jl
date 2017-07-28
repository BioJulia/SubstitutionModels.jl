abstract type K80 <: NASM end


const _πACGT(mod::K80) = 0.25


function P(mod::K80, t::Array{Float64})
  return [P(mod, i) for i in t]
end
