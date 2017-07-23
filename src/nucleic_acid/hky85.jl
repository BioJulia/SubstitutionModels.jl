abstract type HKY85 <: NASM end

struct HKY85abs <: HKY85
  α::Float64
  β::Float64
  πA::Float64
  πC::Float64
  πG::Float64
  πT::Float64
  function HKY85abs(α::Float64, β::Float64,
                  πA::Float64, πC::Float64, πG::Float64, πT::Float64)
    if α <= 0.
      error("HKY85 parameter α must be positive")
    elseif β <= 0.
      error("HKY85 parameter β must be positive")
    elseif sum([πA,πC,πG,πT]) != 1.0
      error("HKY85 frequencies must sum to 1.0")
    elseif any([πA,πC,πG,πT] .<=0.0)
      error("HKY85 frequencies must be positive")
    end
    new(α, β, πA, πC, πG, πT)
  end
end

struct HKY85rel <: HKY85
  κ::Float64
  πA::Float64
  πC::Float64
  πG::Float64
  πT::Float64
  function HKY85rel(κ::Float64,
                  πA::Float64, πC::Float64, πG::Float64, πT::Float64)
    if κ <= 0.
      error("HKY85 parameter κ must be positive")
    elseif sum([πA,πC,πG,πT]) != 1.0
      error("HKY85 frequencies must sum to 1.0")
    elseif any([πA,πC,πG,πT] .<=0.0)
      error("HKY85 frequencies must be positive")
    end
    new(κ, πA, πC, πG, πT)
  end
end

HKY85(α, β, πA, πC, πG, πT) = HKY85abs(α, β, πA, πC, πG, πT)
HKY85(κ, πA, πC, πG, πT) = HKY85rel(κ, πA, πC, πG, πT)

@inline function _μ(mod::HKY85abs)
  return mod.β
end

@inline function _μ(mod::HKY85rel)
  return 1.0
end

@inline function _πA(mod::HKY85)
  return mod.πA
end

@inline function _πC(mod::HKY85)
  return mod.πC
end

@inline function _πG(mod::HKY85)
  return mod.πG
end

@inline function _πT(mod::HKY85)
  return mod.πT
end

"α = r(T/U → C) = r(C → T/U)"
@inline function _α(mod::HKY85abs)
  return mod.α/mod.β
end

"β = r(T/U → A) = r(A → T/U)"
@inline function _β(mod::HKY85abs)
  return 1.0
end

"γ = r(T/U → G) = r(G → T/U)"
@inline function _γ(mod::HKY85abs)
  return 1.0
end

"δ = r(C → A) = r(A → C)"
@inline function _δ(mod::HKY85abs)
  return 1.0
end

"ϵ = r(C → G) = r(G → C)"
@inline function _ϵ(mod::HKY85abs)
  return 1.0
end

"η = r(A → G) = r(G → A)"
@inline function _η(mod::HKY85abs)
  return mod.α/mod.β
end

"α = r(T/U → C) = r(C → T/U)"
@inline function _α(mod::HKY85rel)
  return mod.κ
end

"β = r(T/U → A) = r(A → T/U)"
@inline function _β(mod::HKY85rel)
  return 1.0
end

"γ = r(T/U → G) = r(G → T/U)"
@inline function _γ(mod::HKY85rel)
  return 1.0
end

"δ = r(C → A) = r(A → C)"
@inline function _δ(mod::HKY85rel)
  return 1.0
end

"ϵ = r(C → G) = r(G → C)"
@inline function _ϵ(mod::HKY85rel)
  return 1.0
end

"η = r(A → G) = r(G → A)"
@inline function _η(mod::HKY85rel)
  return mod.κ
end
