abstract type F84 <: NASM end


const _πA(mod::F84) = mod.πA
const _πC(mod::F84) = mod.πC
const _πG(mod::F84) = mod.πG
const _πT(mod::F84) = mod.πT


function P(mod::F84, t::Array{Float64})
  return [P(mod, i) for i in t]
end
