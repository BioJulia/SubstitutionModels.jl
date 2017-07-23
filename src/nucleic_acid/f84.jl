abstract type F84 <: NASM end

struct F84abs <: F84
  α::Float64
  β::Float64
  πA::Float64
  πC::Float64
  πG::Float64
  πT::Float64
  function F84abs(α::Float64, β::Float64,
                  πA::Float64, πC::Float64, πG::Float64, πT::Float64)
    if α <= 0.
      error("F84 parameter α must be positive")
    elseif β <= 0.
      error("F84 parameter β must be positive")
    elseif sum([πA,πC,πG,πT]) != 1.0
      error("F84 frequencies must sum to 1.0")
    elseif any([πA,πC,πG,πT] .<=0.0)
      error("F84 frequencies must be positive")
    end
    new(α, β, πA, πC, πG, πT)
  end
end

struct F84rel <: F84
  κ::Float64
  πA::Float64
  πC::Float64
  πG::Float64
  πT::Float64
  function F84rel(κ::Float64,
                  πA::Float64, πC::Float64, πG::Float64, πT::Float64)
    if κ <= 0.
      error("F84 parameter κ must be positive")
    elseif sum([πA,πC,πG,πT]) != 1.0
      error("F84 frequencies must sum to 1.0")
    elseif any([πA,πC,πG,πT] .<=0.0)
      error("F84 frequencies must be positive")
    end
    new(κ, πA, πC, πG, πT)
  end
end

F84(α, β, πA, πC, πG, πT) = F84abs(α, β, πA, πC, πG, πT)
F84(κ, πA, πC, πG, πT) = F84rel(κ, πA, πC, πG, πT)

@inline function _μ(mod::F84abs)
  return mod.β
end

@inline function _μ(mod::F84rel)
  return 1.0
end

@inline function _πA(mod::F84)
  return mod.πA
end

@inline function _πC(mod::F84)
  return mod.πC
end

@inline function _πG(mod::F84)
  return mod.πG
end

@inline function _πT(mod::F84)
  return mod.πT
end

"α = r(T/U → C) = r(C → T/U)"
@inline function _α(mod::F84abs)
  return mod.α/mod.β
end

"β = r(T/U → A) = r(A → T/U)"
@inline function _β(mod::F84abs)
  return 1.0
end

"γ = r(T/U → G) = r(G → T/U)"
@inline function _γ(mod::F84abs)
  return 1.0
end

"δ = r(C → A) = r(A → C)"
@inline function _δ(mod::F84abs)
  return 1.0
end

"ϵ = r(C → G) = r(G → C)"
@inline function _ϵ(mod::F84abs)
  return 1.0
end

"η = r(A → G) = r(G → A)"
@inline function _η(mod::F84abs)
  return mod.α/mod.β
end

"α = r(T/U → C) = r(C → T/U)"
@inline function _α(mod::F84rel)
  return mod.κ
end

"β = r(T/U → A) = r(A → T/U)"
@inline function _β(mod::F84rel)
  return 1.0
end

"γ = r(T/U → G) = r(G → T/U)"
@inline function _γ(mod::F84rel)
  return 1.0
end

"δ = r(C → A) = r(A → C)"
@inline function _δ(mod::F84rel)
  return 1.0
end

"ϵ = r(C → G) = r(G → C)"
@inline function _ϵ(mod::F84rel)
  return 1.0
end

"η = r(A → G) = r(G → A)"
@inline function _η(mod::F84rel)
  return mod.κ
end
