abstract type TN93 <: NASM end


const _πA(mod::TN93) = mod.πA
const _πC(mod::TN93) = mod.πC
const _πG(mod::TN93) = mod.πG
const _πT(mod::TN93) = mod.πT


function P(mod::TN93, t::Array{Float64})
  return [P(mod, i) for i in t]
end
