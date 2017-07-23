abstract type GTR <: NASM end

struct GTRabs <: GTR
  α::Float64
  β::Float64
  γ::Float64
  δ::Float64
  ϵ::Float64
  η::Float64
  πA::Float64
  πC::Float64
  πG::Float64
  πT::Float64
  function GTRabs(α::Float64, β::Float64, γ::Float64,
                  δ::Float64, ϵ::Float64, η::Float64,
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
    elseif η <= 0.
      error("GTR parameter η must be positive")
    elseif sum([πA,πC,πG,πT]) != 1.0
      error("GTR frequencies must sum to 1.0")
    elseif any([πA,πC,πG,πT] .<= 0.0)
      error("GTR frequencies must be positive")
    end
    new(α, β, γ, δ, ϵ, η, πA, πC, πG, πT)
  end
end

struct GTRrel <: GTR
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

GTR(α, β, γ, δ, ϵ, η, πA, πC, πG, πT) = GTRabs(α, β, γ, δ, ϵ, η, πA, πC, πG, πT)
GTR(α, β, γ, δ, ϵ, πA, πC, πG, πT) = GTRrel(α, β, γ, δ, ϵ, πA, πC, πG, πT)

@inline function _μ(mod::GTR)
  return 1.0
end

@inline function _πA(mod::GTR)
  return mod.πA
end

@inline function _πC(mod::GTR)
  return mod.πC
end

@inline function _πG(mod::GTR)
  return mod.πG
end

@inline function _πT(mod::GTR)
  return mod.πT
end

"α = r(T/U → C) = r(C → T/U)"
@inline function _α(mod::GTR)
  return mod.α
end

"β = r(T/U → A) = r(A → T/U)"
@inline function _β(mod::GTR)
  return mod.β
end

"γ = r(T/U → G) = r(G → T/U)"
@inline function _γ(mod::GTR)
  return mod.β
end

"δ = r(C → A) = r(A → C)"
@inline function _δ(mod::GTR)
  return mod.δ
end

"ϵ = r(C → G) = r(G → C)"
@inline function _ϵ(mod::GTR)
  return mod.ϵ
end

"η = r(A → G) = r(G → A)"
@inline function _η(mod::GTRabs)
  return mod.η
end

"η = r(A → G) = r(G → A)"
@inline function _η(mod::GTRrel)
  return 1.0
end
