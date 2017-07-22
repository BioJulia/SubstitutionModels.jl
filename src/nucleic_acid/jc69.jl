abstract type JC69 <: NASM end

struct JC69abs <: JC69
  λ::Float64
  function JC69abs(λ::Float64)
    if λ <= 0
      error("JC69 parameter λ must be positive")
    end
    new(λ)
  end
end

struct JC69rel <: JC69
end

JC69(λ) = JC69abs(λ)
JC69() = JC69rel()

@inline function μ(mod::JC69abs)
  return mod.λ
end

@inline function μ(mod::JC69rel)
  return 1.0
end

@inline function π(mod::JC69)
  return SVector(0.25, 0.25, 0.25, 0.25)
end

@inline function _α(mod::JC69)
  return 1.
end

@inline function _β(mod::JC69)
  return 1.
end

@inline function _γ(mod::JC69)
  return 1.
end

@inline function _δ(mod::JC69)
  return 1.
end

@inline function _ϵ(mod::JC69)
  return 1.
end

@inline function _ζ(mod::JC69)
  return 1.
end
