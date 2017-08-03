abstract type K80 <: NASM end


function P(mod::K80, t::Array{Float64})
  return [P(mod, i) for i in t]
end
