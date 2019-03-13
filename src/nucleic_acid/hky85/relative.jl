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


HKY85(κ, πA, πC, πG, πT) = HKY85rel(κ, πA, πC, πG, πT)


@inline function Q(mod::HKY85rel)
  κ = mod.κ
  πA = _πA(mod); πC = _πC(mod); πG = _πG(mod); πT = _πT(mod)

  Q₁  = πA
  Q₂  = κ * πA
  Q₃  = πC
  Q₄  = κ * πC
  Q₅  = κ * πG
  Q₆  = πG
  Q₇  = πT
  Q₈  = κ * πT
  Q₉  = -(Q₃ + Q₅ + Q₇)
  Q₁₀ = -(Q₁ + Q₆ + Q₈)
  Q₁₁ = -(Q₂ + Q₃ + Q₇)
  Q₁₂ = -(Q₁ + Q₄ + Q₆)

  return Qmatrix(Q₉,  Q₁,  Q₂,  Q₁,
                 Q₃,  Q₁₀, Q₃,  Q₄,
                 Q₅,  Q₆,  Q₁₁, Q₆,
                 Q₇,  Q₈,  Q₇,  Q₁₂)
end


@inline function P(mod::HKY85rel, t::Float64)
  if t < 0.0
    error("t must be positive")
  end
  κ = mod.κ
  πA = _πA(mod); πC = _πC(mod); πG = _πG(mod); πT = _πT(mod)
  πR = _πR(mod); πY = _πY(mod)

  e₁ = exp(-t)
  e₂ = exp(-(πR * κ + πY) * t)
  e₃ = exp(-(πY * κ + πR) * t)

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
