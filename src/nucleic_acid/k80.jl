abstract type K80 <: NASM end

struct K80abs <: K80
  α::Float64
  β::Float64
  function K80abs(α::Float64, β::Float64)
    if α <= 0.
      error("K80 parameter α must be positive")
    elseif β <= 0.
      error("K80 parameter β must be positive")
    end
    new(α, β)
  end
end

struct K80rel <: K80
  κ::Float64
  function K80rel(κ::Float64)
    if κ <= 0.
      error("K80 parameter κ must be positive")
    end
    new(κ)
  end
end

K80(α, β) = K80abs(α, β)
K80(κ) = K80rel(κ)

@inline function _μ(mod::K80abs)
  return mod.β
end

@inline function _μ(mod::K80rel)
  return 1.0
end

@inline function _πACGT(mod::K80)
  return 0.25
end

"α = r(T/U → C) = r(C → T/U)"
@inline function _α(mod::K80abs)
  return mod.α/mod.β
end

"β = r(T/U → A) = r(A → T/U)"
@inline function _β(mod::K80)
  return 1.0
end

"γ = r(T/U → G) = r(G → T/U)"
@inline function _γ(mod::K80)
  return 1.0
end

"δ = r(C → A) = r(A → C)"
@inline function _δ(mod::K80)
  return 1.0
end

"ϵ = r(C → G) = r(G → C)"
@inline function _ϵ(mod::K80)
  return 1.0
end

"η = r(A → G) = r(G → A)"
@inline function _η(mod::K80abs)
  return mod.α/mod.β
end

"α = r(T/U → C) = r(C → T/U)"
@inline function _α(mod::K80rel)
  return mod.κ
end

"η = r(A → G) = r(G → A)"
@inline function _η(mod::K80rel)
  return mod.κ
end
