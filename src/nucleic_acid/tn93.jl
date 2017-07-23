abstract type TN93 <: NASM end

struct TN93abs <: TN93
  α1::Float64
  α2::Float64
  β::Float64
  πAT::Float64
  πCG::Float64
  function TN93abs(α1::Float64, α2::Float64, β::Float64,
                   πAT::Float64, πCG::Float64)
    if α1 <= 0.
      error("TN93 parameter α1 must be positive")
    elseif α2 <= 0.
      error("TN93 parameter α2 must be positive")
    elseif β <= 0.
      error("TN93 parameter β must be positive")
    elseif sum([πAT,πCG]) != 0.5
      error("TN93 frequencies must sum to 0.5")
    elseif any([πAT,πCG] .<= 0.0)
      error("TN93 frequencies must be positive")
    end
    new(α1, α2, β, πAT, πCG)
  end
end

struct TN93rel <: TN93
  κ1::Float64
  κ2::Float64
  πAT::Float64
  πCG::Float64
  function TN93rel(κ1::Float64, κ2::Float64,
                   πAT::Float64, πCG::Float64)
    if κ1 <= 0.
      error("TN93 parameter κ1 must be positive")
    elseif κ2 <= 0.
      error("TN93 parameter κ2 must be positive")
    elseif sum([πAT,πCG]) != 0.5
      error("TN93 frequencies must sum to 0.5")
    elseif any([πAT,πCG] .<= 0.0)
      error("TN93 frequencies must be positive")
    end
    new(κ1, κ2, πAT, πCG)
  end
end

TN93(α1, α2, β, πAT, πCG) = TN93abs(α1, α2, β, πAT, πCG)
TN93(κ1, κ2, πAT, πCG) = TN93rel(κ1, κ2, πAT, πCG)

@inline function _μ(mod::TN93abs)
  return mod.β
end

@inline function _μ(mod::TN93rel)
  return 1.0
end

@inline function _πAT(mod::TN93)
  return mod.πAT
end

@inline function _πCG(mod::TN93)
  return mod.πCG
end

"α = r(T/U → C) = r(C → T/U)"
@inline function _α(mod::TN93abs)
  return mod.α1/mod.β
end

"β = r(T/U → A) = r(A → T/U)"
@inline function _β(mod::TN93abs)
  return 1.0
end

"γ = r(T/U → G) = r(G → T/U)"
@inline function _γ(mod::TN93abs)
  return 1.0
end

"δ = r(C → A) = r(A → C)"
@inline function _δ(mod::TN93abs)
  return 1.0
end

"ϵ = r(C → G) = r(G → C)"
@inline function _ϵ(mod::TN93abs)
  return 1.0
end

"η = r(A → G) = r(G → A)"
@inline function _η(mod::TN93abs)
  return mod.α2/mod.β
end

"α = r(T/U → C) = r(C → T/U)"
@inline function _α(mod::TN93rel)
  return mod.κ1
end

"β = r(T/U → A) = r(A → T/U)"
@inline function _β(mod::TN93rel)
  return 1.0
end

"γ = r(T/U → G) = r(G → T/U)"
@inline function _γ(mod::TN93rel)
  return 1.0
end

"δ = r(C → A) = r(A → C)"
@inline function _δ(mod::TN93rel)
  return 1.0
end

"ϵ = r(C → G) = r(G → C)"
@inline function _ϵ(mod::TN93rel)
  return 1.0
end

"η = r(A → G) = r(G → A)"
@inline function _η(mod::TN93rel)
  return mod.κ2
end
