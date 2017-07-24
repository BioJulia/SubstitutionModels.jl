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

function show(io::IO, object::K80abs)
  print(io, "\r\e[0m\e[1mK\e[0mimura 19\e[1m80\e[0m model (absolute rate form)
α = $(object.α), β = $(object.β)")
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

function show(io::IO, object::K80rel)
  print(io, "\r\e[0m\e[1mK\e[0mimura 19\e[1m80\e[0m model (relative rate form)
κ = $(object.κ)")
end

K80(α, β) = K80abs(α, β)
K80(κ) = K80rel(κ)

@inline function _μ(mod::K80)
  return 1.0
end

@inline function _πACGT(mod::K80)
  return 0.25
end

"α = r(T/U → C) = r(C → T/U)"
@inline function _α(mod::K80abs)
  return mod.α
end

"β = r(T/U → A) = r(A → T/U)"
@inline function _β(mod::K80abs)
  return mod.β
end

"γ = r(T/U → G) = r(G → T/U)"
@inline function _γ(mod::K80abs)
  return mod.β
end

"δ = r(C → A) = r(A → C)"
@inline function _δ(mod::K80abs)
  return mod.β
end

"ϵ = r(C → G) = r(G → C)"
@inline function _ϵ(mod::K80abs)
  return mod.β
end

"η = r(A → G) = r(G → A)"
@inline function _η(mod::K80abs)
  return mod.α
end

"α = r(T/U → C) = r(C → T/U)"
@inline function _α(mod::K80rel)
  return mod.κ
end

"β = r(T/U → A) = r(A → T/U)"
@inline function _β(mod::K80rel)
  return 1.0
end

"γ = r(T/U → G) = r(G → T/U)"
@inline function _γ(mod::K80rel)
  return 1.0
end

"δ = r(C → A) = r(A → C)"
@inline function _δ(mod::K80rel)
  return 1.0
end

"ϵ = r(C → G) = r(G → C)"
@inline function _ϵ(mod::K80rel)
  return 1.0
end

"η = r(A → G) = r(G → A)"
@inline function _η(mod::K80rel)
  return mod.κ
end


"Generate a P matrix for a `K80abs`, of the form:

  [[A→A, A→C, A→G, A→T]
   [C→A, C→C, C→G, C→T]
   [G→A, G→C, G→G, G→T]
   [T→A, T→C, T→G, T→T]]

for a specified time"
function P(mod::K80abs, t::Float64)
  α = mod.α
  β = mod.β
  ω = 0.25 * exp(-4 * β * t)
  τ = 0.5 * exp(-2 * (α + β) * t)
  P₁ = 0.25 + ω + τ
  P₂ = 0.25 + ω - τ
  P₃ = 0.25 - ω
  return SMatrix{4, 4, Float64}(P₁, P₃, P₂, P₃,
                                P₃, P₁, P₃, P₂,
                                P₂, P₃, P₁, P₃,
                                P₃, P₂, P₃, P₁)
end

"Generate a P matrix for a `K80rel`, of the form:

  [[A→A, A→C, A→G, A→T]
   [C→A, C→C, C→G, C→T]
   [G→A, G→C, G→G, G→T]
   [T→A, T→C, T→G, T→T]]

for a specified time"
function P(mod::K80rel, t::Float64)
  κ = mod.κ
  ω = 0.25 * exp(-4 * t/(κ + 2))
  τ = 0.5 * exp(-2 * t * (κ + 1)/(κ + 2))
  P₁ = 0.25 + ω + τ
  P₂ = 0.25 + ω - τ
  P₃ = 0.25 - ω
  return SMatrix{4, 4, Float64}(P₁, P₃, P₂, P₃,
                                P₃, P₁, P₃, P₂,
                                P₂, P₃, P₁, P₃,
                                P₃, P₂, P₃, P₁)
end

"Generate an array of P matrices for a specified array of times"
function P(mod::K80, t::Array{Float64})
  return [P(mod, i) for i in t]
end
