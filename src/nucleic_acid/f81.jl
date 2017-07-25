abstract type F81 <: NASM end

struct F81abs <: F81
  β::Float64
  πA::Float64
  πC::Float64
  πG::Float64
  πT::Float64
  function F81abs(β::Float64,
                  πA::Float64, πC::Float64, πG::Float64, πT::Float64)
    if β <= 0.
      error("F81 parameter β must be positive")
    elseif sum([πA,πC,πG,πT]) != 1.0
      error("F81 frequencies must sum to 1.0")
    elseif any([πA,πC,πG,πT] .<=0.0)
      error("F81 frequencies must be positive")
    end
    new(β, πA, πC, πG, πT)
  end
end

function show(io::IO, object::F81abs)
  print(io, "\r\e[0m\e[1mF\e[0melsenstein 19\e[1m81\e[0m model (absolute rate form)
β = $(object.β), π = [$(object.πA), $(object.πC), $(object.πG), $(object.πT)]")
end

struct F81rel <: F81
  πA::Float64
  πC::Float64
  πG::Float64
  πT::Float64
  function F81rel(πA::Float64, πC::Float64, πG::Float64, πT::Float64)
    if sum([πA,πC,πG,πT]) != 1.0
      error("F81 frequencies must sum to 1.0")
    elseif any([πA,πC,πG,πT] .<=0.0)
      error("F81 frequencies must be positive")
    end
    new(πA, πC, πG, πT)
  end
end

function show(io::IO, object::F81rel)
  print(io, "\r\e[0m\e[1mF\e[0melsenstein 19\e[1m81\e[0m model (relative rate form)
π = [$(object.πA), $(object.πC), $(object.πG), $(object.πT)]")
end

F81(β, πA, πC, πG, πT) = F81abs(β, πA, πC, πG, πT)
F81(πA, πC, πG, πT) = F81rel(πA, πC, πG, πT)

@inline function _μ(mod::F81abs)
  return mod.β
end

@inline function _μ(mod::F81rel)
  return 1.0
end

@inline function _πA(mod::F81)
  return mod.πA
end

@inline function _πC(mod::F81)
  return mod.πC
end

@inline function _πG(mod::F81)
  return mod.πG
end

@inline function _πT(mod::F81)
  return mod.πT
end

"α = r(T/U → C) = r(C → T/U)"
@inline function _α(mod::F81)
  return 1.0
end

"β = r(T/U → A) = r(A → T/U)"
@inline function _β(mod::F81)
  return 1.0
end

"γ = r(T/U → G) = r(G → T/U)"
@inline function _γ(mod::F81)
  return 1.0
end

"δ = r(C → A) = r(A → C)"
@inline function _δ(mod::F81)
  return 1.0
end

"ϵ = r(C → G) = r(G → C)"
@inline function _ϵ(mod::F81)
  return 1.0
end

"η = r(A → G) = r(G → A)"
@inline function _η(mod::F81)
  return 1.0
end

"Generate a P matrix for a `F84`, of the form:

  [[A→A, A→C, A→G, A→T]
   [C→A, C→C, C→G, C→T]
   [G→A, G→C, G→G, G→T]
   [T→A, T→C, T→G, T→T]]

for a specified time"
@inline function P(mod::F81, t::Float64)
  πA = _πA(mod)
  πC = _πC(mod)
  πG = _πG(mod)
  πT = _πT(mod)

  πR = _πR(mod)
  πY = _πY(mod)

  β = _μ(mod)
  ω = exp(-β * t)

  P₁  = πA + ((πA * πY / πR) + (πG / πR)) * ω
  P₂  = πC + ((πC * πR / πY) + (πT / πY)) * ω
  P₃  = πG + ((πG * πY / πR) + (πA / πR)) * ω
  P₄  = πT + ((πT * πR / πY) + (πC / πY)) * ω
  P₅  = πA * (1.0 - ω)
  P₆  = πC * (1.0 - ω)
  P₇  = πG * (1.0 - ω)
  P₈  = πT * (1.0 - ω)

  return SMatrix{4, 4, Float64}(P₁,  P₅,  P₅,  P₅,
                                P₆,  P₂,  P₆,  P₆,
                                P₇,  P₇,  P₃,  P₇,
                                P₈,  P₈,  P₈,  P₄)
end

"Generate an array of P matrices for a specified array of times"
function P(mod::F81, t::Array{Float64})
  return [P(mod, i) for i in t]
end
