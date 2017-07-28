abstract type TN93 <: NASM end


const _πAT(mod::TN93) = mod.πAT
const _πCG(mod::TN93) = mod.πCG


function P(mod::TN93, t::Array{Float64})
  return [P(mod, i) for i in t]
end
