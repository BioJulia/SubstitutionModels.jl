abstract type F84 <: NASM end

struct F84abs <: F84
  κ::Float64
  β::Float64
  πA::Float64
  πC::Float64
  πG::Float64
  πT::Float64
  function F84abs(κ::Float64, β::Float64,
                  πA::Float64, πC::Float64, πG::Float64, πT::Float64)
    if κ <= 0.
      error("F84 parameter α must be positive")
    elseif β <= 0.
      error("F84 parameter β must be positive")
    elseif sum([πA,πC,πG,πT]) != 1.0
      error("F84 frequencies must sum to 1.0")
    elseif any([πA,πC,πG,πT] .<= 0.0)
      error("F84 frequencies must be positive")
    end
    new(κ, β, πA, πC, πG, πT)
  end
end

function show(io::IO, object::F84abs)
  print(io, "\r\e[0m\e[1mF\e[0melsenstein 19\e[1m84\e[0m substitution model (absolute rate form)
κ = $(object.κ), β = $(object.β), π = [$(object.πA), $(object.πC), $(object.πG), $(object.πT)]")
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

function show(io::IO, object::F84rel)
  print(io, "\r\e[0m\e[1mF\e[0melsenstein 19\e[1m84\e[0m substitution model (relative rate form)
κ = $(object.κ), π = [$(object.πA), $(object.πC), $(object.πG), $(object.πT)]")
end

F84(κ, β, πA, πC, πG, πT) = F84abs(κ, β, πA, πC, πG, πT)
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
@inline function _α(mod::F84)
  return 1.0 + mod.κ/_πR(mod)
end

"β = r(T/U → A) = r(A → T/U)"
@inline function _β(mod::F84)
  return 1.0
end

"γ = r(T/U → G) = r(G → T/U)"
@inline function _γ(mod::F84)
  return 1.0
end

"δ = r(C → A) = r(A → C)"
@inline function _δ(mod::F84)
  return 1.0
end

"ϵ = r(C → G) = r(G → C)"
@inline function _ϵ(mod::F84)
  return 1.0
end

"η = r(A → G) = r(G → A)"
@inline function _η(mod::F84)
  return 1.0 + mod.κ/_πY(mod)
end

"Generate a P matrix for a `F84`, of the form:

  [[A→A, A→C, A→G, A→T]
   [C→A, C→C, C→G, C→T]
   [G→A, G→C, G→G, G→T]
   [T→A, T→C, T→G, T→T]]

for a specified time"
@inline function P(mod::F84, t::Float64)
  πA = _πA(mod)
  πC = _πC(mod)
  πG = _πG(mod)
  πT = _πT(mod)

  πR = _πR(mod)
  πY = _πY(mod)

  β = _μ(mod)
  κ = mod.κ

  α₁ = (1.0 + κ / πY) * β
  α₂ = (1.0 + κ / πR) * β

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

"Generate an array of P matrices for a specified array of times"
function P(mod::F84, t::Array{Float64})
  return [P(mod, i) for i in t]
end
