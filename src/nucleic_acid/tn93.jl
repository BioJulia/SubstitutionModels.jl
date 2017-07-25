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

function show(io::IO, object::TN93abs)
  print(io, "\r\e[0m\e[1mT\e[0mamura and \e[1mN\e[0mei 19\e[1m93\e[0m model (absolute rate form)
α1 = $(object.α1), α2 = $(object.α1), β = $(object.β), π = [$(object.πAT), $(object.πCG), $(object.πCG), $(object.πAT)]")
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

function show(io::IO, object::TN93rel)
  print(io, "\r\e[0m\e[1mT\e[0mamura and \e[1mN\e[0mei 19\e[1m93\e[0m model (relative rate form)
κ1 = $(object.κ1), κ2 = $(object.κ1), π = [$(object.πAT), $(object.πCG), $(object.πCG), $(object.πAT)]")
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


"Generate a P matrix for a `TN93abs`, of the form:

  [[A→A, A→C, A→G, A→T]
   [C→A, C→C, C→G, C→T]
   [G→A, G→C, G→G, G→T]
   [T→A, T→C, T→G, T→T]]

for a specified time"
@inline function P(mod::TN93abs, t::Float64)
  if t < 0.0
    error("t must be positive")
  end
  πA = _πA(mod)
  πC = _πC(mod)
  πG = _πG(mod)
  πT = _πT(mod)

  πR = _πR(mod)
  πY = _πY(mod)

  β = mod.β
  α₁ = mod.α1
  α₂ = mod.α2

  e₁ = exp(-β * t)
  e₂ = exp(-(πR * α₂ + πY * β) * t)
  e₃ = exp(-(πY * α₁ + πR * β) * t)

  P₁  = πA + (πA * πY / πR) * e₁ + (πG / πR) * e₂
  P₂  = πC + (πT * πR / πY) * e₁ + (πT / πY) * e₃
  P₃  = πG + (πG * πY / πR) * e₁ + (πA / πR) * e₂
  P₄  = πT + (πT * πR / πY) * e₁ + (πC / πY) * e₃
  P₅  = πA * (1 - e₁)
  P₆  = πA + (πA * πY / πR) * e₁ - (πA / πR) * e₂
  P₇  = πC * (1 - e₁)
  P₈  = πC + (πT * πR / πY) * e₁ - (πC / πY) * e₃
  P₉  = πG + (πG * πY / πR) * e₁ - (πG / πR) * e₂
  P₁₀ = πG * (1 - e₁)
  P₁₁ = πT * (1 - e₁)
  P₁₂ = πT + (πT * πR / πY) * e₁ - (πT / πY) * e₃

  return SMatrix{4, 4, Float64}(P₁,  P₅,  P₆,  P₅,
                                P₇,  P₂,  P₇,  P₈,
                                P₉,  P₁₀, P₃,  P₁₀,
                                P₁₁, P₁₂, P₁₁, P₄)
end

"Generate a P matrix for a `TN93rel`, of the form:

  [[A→A, A→C, A→G, A→T]
   [C→A, C→C, C→G, C→T]
   [G→A, G→C, G→G, G→T]
   [T→A, T→C, T→G, T→T]]

for a specified time"
@inline function P(mod::TN93rel, t::Float64)
  if t < 0.0
    error("t must be positive")
  end
  πA = _πA(mod)
  πC = _πC(mod)
  πG = _πG(mod)
  πT = _πT(mod)

  πR = _πR(mod)
  πY = _πY(mod)

  κ₁ = mod.κ1
  κ₂ = mod.κ2

  e₁ = exp(-t)
  e₂ = exp(-(πR * κ₂ + πY) * t)
  e₃ = exp(-(πY * κ₁ + πR) * t)

  P₁  = πA + (πA * πY / πR) * e₁ + (πG / πR) * e₂
  P₂  = πC + (πT * πR / πY) * e₁ + (πT / πY) * e₃
  P₃  = πG + (πG * πY / πR) * e₁ + (πA / πR) * e₂
  P₄  = πT + (πT * πR / πY) * e₁ + (πC / πY) * e₃
  P₅  = πA * (1 - e₁)
  P₆  = πA + (πA * πY / πR) * e₁ - (πA / πR) * e₂
  P₇  = πC * (1 - e₁)
  P₈  = πC + (πT * πR / πY) * e₁ - (πC / πY) * e₃
  P₉  = πG + (πG * πY / πR) * e₁ - (πG / πR) * e₂
  P₁₀ = πG * (1 - e₁)
  P₁₁ = πT * (1 - e₁)
  P₁₂ = πT + (πT * πR / πY) * e₁ - (πT / πY) * e₃

  return SMatrix{4, 4, Float64}(P₁,  P₅,  P₆,  P₅,
                                P₇,  P₂,  P₇,  P₈,
                                P₉,  P₁₀, P₃,  P₁₀,
                                P₁₁, P₁₂, P₁₁, P₄)
end

"Generate an array of P matrices for a specified array of times"
function P(mod::TN93, t::Array{Float64})
  return [P(mod, i) for i in t]
end
