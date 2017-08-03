abstract type F81 <: NASM end

const _πA(mod::F81) = mod.πA
const _πC(mod::F81) = mod.πC
const _πG(mod::F81) = mod.πG
const _πT(mod::F81) = mod.πT


function P(mod::F81, t::Array{Float64})
  return [P(mod, i) for i in t]
end
