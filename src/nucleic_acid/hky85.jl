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
    elseif any([πA,πC,πG,πT] .<= 0.0)
      error("HKY85 frequencies must be positive")
    end
    new(α, β, πA, πC, πG, πT)
  end
end

function show(io::IO, object::HKY85abs)
  print(io, "\r\e[0m\e[1mH\e[0masegawa, \e[1mK\e[0mishino, and \e[1mY\e[0mano 19\e[1m85\e[0m model (absolute rate form)
α = $(object.α), β = $(object.β), π = [$(object.πA), $(object.πC), $(object.πG), $(object.πT)]")
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

function show(io::IO, object::HKY85rel)
  print(io, "\r\e[0m\e[1mH\e[0masegawa, \e[1mK\e[0mishino, and \e[1mY\e[0mano 19\e[1m85\e[0m model (absolute rate form)
κ = $(object.κ), π = [$(object.πA), $(object.πC), $(object.πG), $(object.πT)]")
end

HKY85(α, β, πA, πC, πG, πT) = HKY85abs(α, β, πA, πC, πG, πT)
HKY85(κ, πA, πC, πG, πT) = HKY85rel(κ, πA, πC, πG, πT)

@inline function _μ(mod::HKY85)
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
  return mod.α
end

"β = r(T/U → A) = r(A → T/U)"
@inline function _β(mod::HKY85abs)
  return mod.β
end

"γ = r(T/U → G) = r(G → T/U)"
@inline function _γ(mod::HKY85abs)
  return mod.β
end

"δ = r(C → A) = r(A → C)"
@inline function _δ(mod::HKY85abs)
  return mod.β
end

"ϵ = r(C → G) = r(G → C)"
@inline function _ϵ(mod::HKY85abs)
  return mod.β
end

"η = r(A → G) = r(G → A)"
@inline function _η(mod::HKY85abs)
  return mod.α
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

"Generate a P matrix for a `HKY85abs`, of the form:

  [[A→A, A→C, A→G, A→T]
   [C→A, C→C, C→G, C→T]
   [G→A, G→C, G→G, G→T]
   [T→A, T→C, T→G, T→T]]

for a specified time"
@inline function P(mod::HKY85abs, t::Float64)
  πA = _πA(mod)
  πC = _πC(mod)
  πG = _πG(mod)
  πT = _πT(mod)

  πR = _πR(mod)
  πY = _πY(mod)

  α = mod.α
  β = mod.β

  e₁ = exp(-β * t)
  e₂ = exp(-(πR * α + πY * β) * t)
  e₃ = exp(-(πY * α + πR * β) * t)

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

"Generate a P matrix for a `HKY85rel`, of the form:

  [[A→A, A→C, A→G, A→T]
   [C→A, C→C, C→G, C→T]
   [G→A, G→C, G→G, G→T]
   [T→A, T→C, T→G, T→T]]

for a specified time"
@inline function P(mod::HKY85rel, t::Float64)
  πA = _πA(mod)
  πC = _πC(mod)
  πG = _πG(mod)
  πT = _πT(mod)

  πR = _πR(mod)
  πY = _πY(mod)

  κ = mod.κ

  e₁ = exp(-t)
  e₂ = exp(-(πR * κ + πY) * t)
  e₃ = exp(-(πY * κ + πR) * t)

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
function P(mod::HKY85, t::Array{Float64})
  return [P(mod, i) for i in t]
end
