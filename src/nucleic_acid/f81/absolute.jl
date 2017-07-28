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


F81(β, πA, πC, πG, πT) = F81abs(β, πA, πC, πG, πT)


@inline function Q(mod::F81abs)
  β = mod.β
  πA = _πA(mod); πC = _πC(mod); πG = _πG(mod); πT = _πT(mod)

  Q₁ = β * πA
  Q₂ = β * πC
  Q₃ = β * πG
  Q₄ = β * πT
  Q₅ = -(Q₂ + Q₃ + Q₄)
  Q₆ = -(Q₁ + Q₃ + Q₄)
  Q₇ = -(Q₁ + Q₂ + Q₄)
  Q₈ = -(Q₁ + Q₂ + Q₃)

  return Qmatrix(Q₅, Q₁, Q₁, Q₁,
                 Q₂, Q₆, Q₂, Q₂,
                 Q₃, Q₃, Q₇, Q₃,
                 Q₄, Q₄, Q₄, Q₈)
end


@inline function P(mod::F81abs, t::Float64)
  if t < 0.0
    error("t must be positive")
  end
  β = mod.β
  πA = _πA(mod); πC = _πC(mod); πG = _πG(mod); πT = _πT(mod)
  πR = πA + πG; πY = πT + πC

  e₁ = exp(-β * t)

  P₁  = πA + ((πA * πY / πR) + (πG / πR)) * e₁
  P₂  = πC + ((πC * πR / πY) + (πT / πY)) * e₁
  P₃  = πG + ((πG * πY / πR) + (πA / πR)) * e₁
  P₄  = πT + ((πT * πR / πY) + (πC / πY)) * e₁
  P₅  = πA * (1.0 - e₁)
  P₆  = πC * (1.0 - e₁)
  P₇  = πG * (1.0 - e₁)
  P₈  = πT * (1.0 - e₁)

  return Pmatrix(P₁,  P₅,  P₅,  P₅,
                 P₆,  P₂,  P₆,  P₆,
                 P₇,  P₇,  P₃,  P₇,
                 P₈,  P₈,  P₈,  P₄)
end
