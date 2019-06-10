struct F84abs <: F84
  κ::Float64
  β::Float64
  πA::Float64
  πC::Float64
  πG::Float64
  πT::Float64
  function F84abs(κ::Float64, β::Float64,
                  πA::Float64, πC::Float64, πG::Float64, πT::Float64,
                  safe::Bool=true)
    if safe
      if κ <= 0.
        error("F84 parameter α must be positive")
      elseif β <= 0.
        error("F84 parameter β must be positive")
      elseif sum([πA,πC,πG,πT]) != 1.0
        error("F84 frequencies must sum to 1.0")
      elseif any([πA,πC,πG,πT] .<= 0.0)
        error("F84 frequencies must be positive")
      end
    end
    new(κ, β, πA, πC, πG, πT)
  end
end


function Base.show(io::IO, object::F84abs)
  print(io, "\r\e[0m\e[1mF\e[0melsenstein 19\e[1m84\e[0m substitution model (absolute rate form)
κ = $(object.κ), β = $(object.β), π = [$(object.πA), $(object.πC), $(object.πG), $(object.πT)]")
end


@inline function Q(mod::F84abs)
  β = mod.β; κ = mod.κ
  πA = _πA(mod); πC = _πC(mod); πG = _πG(mod); πT = _πT(mod)
  πR = _πR(mod); πY = _πY(mod)

  α₁ = (1.0 + κ / πY) * β
  α₂ = (1.0 + κ / πR) * β

  Q₁  = β * πA
  Q₂  = α₂ * πA
  Q₃  = β * πC
  Q₄  = α₁ * πC
  Q₅  = α₂ * πG
  Q₆  = β * πG
  Q₇  = β * πT
  Q₈  = α₁ * πT
  Q₉  = -(Q₃ + Q₅ + Q₇)
  Q₁₀ = -(Q₁ + Q₆ + Q₈)
  Q₁₁ = -(Q₂ + Q₃ + Q₇)
  Q₁₂ = -(Q₁ + Q₄ + Q₆)

  return Qmatrix(Q₉,  Q₁,  Q₂,  Q₁,
                 Q₃,  Q₁₀, Q₃,  Q₄,
                 Q₅,  Q₆,  Q₁₁, Q₆,
                 Q₇,  Q₈,  Q₇,  Q₁₂)
end


@inline function P(mod::F84abs, t::Float64)
  if t < 0.0
    error("t must be positive")
  end
  β = mod.β; κ = mod.κ
  πA = _πA(mod); πC = _πC(mod); πG = _πG(mod); πT = _πT(mod)
  πR = _πR(mod); πY = _πY(mod)

  α₁ = (1.0 + κ / πY) * β
  α₂ = (1.0 + κ / πR) * β

  e₁ = exp(-β * t)
  e₂ = exp(-(πR * α₂ + πY * β) * t)
  e₃ = exp(-(πY * α₁ + πR * β) * t)

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
