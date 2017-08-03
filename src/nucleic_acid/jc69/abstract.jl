abstract type JC69 <: NASM end


function P(mod::JC69, t::Array{Float64})
  return [P(mod, i) for i in t]
end
