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

function show(io::IO, object::JC69abs)
  print(io, "\r\e[0m\e[1mJ\e[0mukes and \e[1mC\e[0mantor 19\e[1m69\e[0m model (absolute rate form)
λ = $(object.λ)")
end

struct JC69rel <: JC69
end

function show(io::IO, object::JC69rel)
  print(io, "\r\e[0m\e[1mJ\e[0mukes and \e[1mC\e[0mantor 19\e[1m69\e[0m model (relative rate form)")
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

"Generate a P matrix for a `JC69` model, of the form:

  [[A→A, A→C, A→G, A→T]
   [C→A, C→C, C→G, C→T]
   [G→A, G→C, G→G, G→T]
   [T→A, T→C, T→G, T→T]]

for a specified time"
@inline function P(mod::JC69, t::Float64)
  if t < 0.0
    error("t must be positive")
  end
  μ = _μ(mod)
  ω = exp(-t * μ)
  P₁ = 0.25 + 0.75 * ω
  P₂ = 0.25 + 0.25 * ω
  return SMatrix{4, 4, Float64}(P₁, P₂, P₂, P₂,
                                P₂, P₁, P₂, P₂,
                                P₂, P₂, P₁, P₂,
                                P₂, P₂, P₂, P₁)
end

"Generate an array of P matrices for a specified array of times"
function P(mod::JC69, t::Array{Float64})
  return [P(mod, i) for i in t]
end
