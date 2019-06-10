struct TN93rel <: TN93
  κ1::Float64
  κ2::Float64
  πA::Float64
  πC::Float64
  πG::Float64
  πT::Float64
  function TN93rel(κ1::Float64, κ2::Float64,
                   πA::Float64, πC::Float64, πG::Float64, πT::Float64,
                   safe::Bool=true)
    if safe
      if κ1 <= 0.
        error("TN93 parameter κ1 must be positive")
      elseif κ2 <= 0.
        error("TN93 parameter κ2 must be positive")
      elseif sum([πA, πC, πG, πT]) != 1.0
        error("TN93 frequencies must sum to 1.0")
      elseif any([πA, πC, πG, πT] .<= 0.0)
        error("TN93 frequencies must be positive")
      end
    end
    new(κ1, κ2, πA, πC, πG, πT)
  end
end


function Base.show(io::IO, object::TN93rel)
  print(io, "\r\e[0m\e[1mT\e[0mamura and \e[1mN\e[0mei 19\e[1m93\e[0m model (relative rate form)
κ1 = $(object.κ1), κ2 = $(object.κ1), π = [$(object.πA), $(object.πC), $(object.πG), $(object.πT)]")
end


@inline function Q(mod::TN93rel)
  κ₁ = mod.κ1; κ₂ = mod.κ2
  πA = _πA(mod); πC = _πC(mod); πG = _πG(mod); πT = _πT(mod)
  πR = _πR(mod); πY = _πY(mod)

  Q₁  = πA
  Q₂  = κ₂ * πA
  Q₃  = πC
  Q₄  = κ₁ * πC
  Q₅  = κ₂ * πG
  Q₆  = πG
  Q₇  = πT
  Q₈  = κ₁ * πT
  Q₉  = -(Q₃ + Q₅ + Q₇)
  Q₁₀ = -(Q₁ + Q₆ + Q₈)
  Q₁₁ = -(Q₂ + Q₃ + Q₇)
  Q₁₂ = -(Q₁ + Q₄ + Q₆)

  return Qmatrix(Q₉,  Q₁,  Q₂,  Q₁,
                 Q₃,  Q₁₀, Q₃,  Q₄,
                 Q₅,  Q₆,  Q₁₁, Q₆,
                 Q₇,  Q₈,  Q₇,  Q₁₂)
end


@inline function P(mod::TN93rel, t::Float64)
  if t < 0.0
    error("t must be positive")
  end
  κ₁ = mod.κ1; κ₂ = mod.κ2
  πA = _πA(mod); πC = _πC(mod); πG = _πG(mod); πT = _πT(mod)
  πR = _πR(mod); πY = _πY(mod)
  e₁ = exp(-t)
  e₂ = exp(-(πR * κ₂ + πY) * t)
  e₃ = exp(-(πY * κ₁ + πR) * t)

  P₁  = πA + (πA * πY / πR) * e₁ + (πG / πR) * e₂
  P₂  = πC + (πC * πR / πY) * e₁ + (πT / πY) * e₃
  P₃  = πG + (πG * πY / πR) * e₁ + (πA / πR) * e₂
  P₄  = πT + (πT * πR / πY) * e₁ + (πC / πY) * e₃
  P₅  = πA * (1 - e₁)
  P₆  = πA + (πA * πY / πR) * e₁ - (πA / πR) * e₂
  P₇  = πC * (1 - e₁)
  P₈  = πC + (πC * πR / πY) * e₁ - (πC / πY) * e₃
  P₉  = πG + (πG * πY / πR) * e₁ - (πG / πR) * e₂
  P₁₀ = πG * (1 - e₁)
  P₁₁ = πT * (1 - e₁)
  P₁₂ = πT + (πT * πR / πY) * e₁ - (πT / πY) * e₃

  return Pmatrix(P₁,  P₅,  P₆,  P₅,
                 P₇,  P₂,  P₇,  P₈,
                 P₉,  P₁₀, P₃,  P₁₀,
                 P₁₁, P₁₂, P₁₁, P₄)
end
