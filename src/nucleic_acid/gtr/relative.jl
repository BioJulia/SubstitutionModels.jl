mutable struct GTRrel <: GTR
  α::Float64
  β::Float64
  γ::Float64
  δ::Float64
  ϵ::Float64
  πA::Float64
  πC::Float64
  πG::Float64
  πT::Float64
  function GTRrel(α::Float64, β::Float64, γ::Float64,
                  δ::Float64, ϵ::Float64,
                  πA::Float64, πC::Float64, πG::Float64, πT::Float64)
    if α <= 0.
      error("GTR parameter α must be positive")
    elseif β <= 0.
      error("GTR parameter β must be positive")
    elseif γ <= 0.
      error("GTR parameter γ must be positive")
    elseif δ <= 0.
      error("GTR parameter δ must be positive")
    elseif ϵ <= 0.
      error("GTR parameter ϵ must be positive")
    elseif sum([πA,πC,πG,πT]) != 1.0
      error("GTR frequencies must sum to 1.0")
    elseif any([πA,πC,πG,πT] .<= 0.0)
      error("GTR frequencies must be positive")
    end
    new(α, β, γ, δ, ϵ, πA, πC, πG, πT)
  end
end

"""
```julia-repl
julia> model = GTR(0.3, 0.6, 0.6, 0.4, 0.4, 0.25, 0.25, 0.25, 0.25);

julia> setrate!(model, [0.5, 0.5, 0.5, 0.5, 0.5, 0.25, 0.25, 0.25, 0.25])
Generalised Time Reversible model (relative rate form)
α = 0.5, β = 0.5, γ = 0.5, δ = 0.5, ϵ = 0.5, π = [0.25, 0.25, 0.25, 0.25]
```
"""
@inline function setrate!(mod::GTRrel, rate::Vector{Float64})
  #check length of vector
  if length(rate) != 9
    @error "GTR rate must be a vector of length 9"
  else 
      #separate vector into pieces
    α = rate[1]
    β = rate[2]
    γ = rate[3]
    δ = rate[4]
    ϵ = rate[5]
    πA = rate[6]
    πC = rate[7]
    πG = rate[8]
    πT = rate[9]
    #check correctness of rates (from above)
    if α <= 0.
      @error "GTR parameter α must be positive"
    elseif β <= 0.
      @error "GTR parameter β must be positive"
    elseif γ <= 0.
      @error "GTR parameter γ must be positive"
    elseif δ <= 0.
      @error "GTR parameter δ must be positive"
    elseif ϵ <= 0.
      @error "GTR parameter ϵ must be positive"
    elseif sum([πA,πC,πG,πT]) != 1.0
      @error "GTR frequencies must sum to 1.0"
    elseif any([πA,πC,πG,πT] .<= 0.0)
      @error "GTR frequencies must be positive"
    end
    #add rate
    mod.α = α
    mod.β = β
    mod.γ = γ
    mod.δ = δ
    mod.ϵ = ϵ
    mod.πA = πA
    mod.πC = πC
    mod.πG = πG
    mod.πT = πT
  end
  return mod
end

function show(io::IO, object::GTRrel)
  print(io, "\r\e[0m\e[1mG\e[0meneralised \e[1mT\e[0mime \e[1mR\e[0meversible model (relative rate form)
α = $(object.α), β = $(object.β), γ = $(object.γ), δ = $(object.δ), ϵ = $(object.ϵ), π = [$(object.πA), $(object.πC), $(object.πG), $(object.πT)]")
end


GTR(α, β, γ, δ, ϵ, πA, πC, πG, πT) = GTRrel(α, β, γ, δ, ϵ, πA, πC, πG, πT)


const _α(mod::GTRrel) = mod.α
const _β(mod::GTRrel) = mod.β
const _γ(mod::GTRrel) = mod.γ
const _δ(mod::GTRrel) = mod.δ
const _ϵ(mod::GTRrel) = mod.ϵ
const _η(mod::GTRrel) = 1.0
