abstract type JC69 <: NASM end

struct JC69abs <: JC69
  λ::Float64
  function JC69abs(λ::Float64)
    if λ <= 0.
      error("JC69 parameter λ must be positive")
    end
    new(λ)
  end
end

struct JC69rel <: JC69
end

JC69(λ) = JC69abs(λ)
JC69() = JC69rel()

@inline function _μ(mod::JC69abs)
  return mod.λ
end

@inline function _μ(mod::JC69rel)
  return 1.0
end

@inline function _πACGT(mod::JC69)
  return 0.25
end

"α = r(T/U → C) = r(C → T/U)"
@inline function _α(mod::JC69)
  return 1.
end

"β = r(T/U → A) = r(A → T/U)"
@inline function _β(mod::JC69)
  return 1.
end

"γ = r(T/U → G) = r(G → T/U)"
@inline function _γ(mod::JC69)
  return 1.
end

"δ = r(C → A) = r(A → C)"
@inline function _δ(mod::JC69)
  return 1.
end

"ϵ = r(C → G) = r(G → C)"
@inline function _ϵ(mod::JC69)
  return 1.
end

"η = r(A → G) = r(G → A)"
@inline function _η(mod::JC69)
  return 1.
end
