mutable struct F84rel <: F84
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

"""
```julia-repl
julia> model = F84(0.6, 0.25, 0.25, 0.25, 0.25);

julia> setrate!(model, [0.5, 0.25, 0.25, 0.25, 0.25])
Felsenstein 1984 substitution model (relative rate form)
κ = 0.5, π = [0.25, 0.25, 0.25, 0.25]
```
"""
@inline function setrate!(mod::F84rel, rate::Vector{Float64})
  #check length of vector
  if length(rate) != 5
    @error "F84 rate must be a vector of length 5"
  else 
    #separate vector into pieces
    κ = rate[1]
    πA = rate[2]
    πC = rate[3]
    πG = rate[4]
    πT = rate[5]
    #check correctness of rates (from above)
    if κ <= 0.
      @error "F84 parameter κ must be positive"
    elseif sum([πA,πC,πG,πT]) != 1.0
      @error "F84 frequencies must sum to 1.0"
    elseif any([πA,πC,πG,πT] .<=0.0)
      @error "F84 frequencies must be positive"
    end
    #assign rates
    mod.κ = κ
    mod.πA = πA
    mod.πC = πC
    mod.πG = πG
    mod.πT = πT
  end
  return mod
end


function show(io::IO, object::F84rel)
  print(io, "\r\e[0m\e[1mF\e[0melsenstein 19\e[1m84\e[0m substitution model (relative rate form)
κ = $(object.κ), π = [$(object.πA), $(object.πC), $(object.πG), $(object.πT)]")
end


F84(κ, πA, πC, πG, πT) = F84rel(κ, πA, πC, πG, πT)


@inline function Q(mod::F84rel)
  κ = mod.κ
  πA = _πA(mod); πC = _πC(mod); πG = _πG(mod); πT = _πT(mod)
  πR = _πR(mod); πY = _πY(mod)

  α₁ = (1.0 + κ / πY)
  α₂ = (1.0 + κ / πR)

  Q₁  = πA
  Q₂  = α₂ * πA
  Q₃  = πC
  Q₄  = α₁ * πC
  Q₅  = α₂ * πG
  Q₆  = πG
  Q₇  = πT
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


@inline function P(mod::F84rel, t::Float64)
  if t < 0.0
    error("t must be positive")
  end
  κ = mod.κ
  πA = _πA(mod); πC = _πC(mod); πG = _πG(mod); πT = _πT(mod)
  πR = _πR(mod); πY = _πY(mod)

  α₁ = (1.0 + κ / πY)
  α₂ = (1.0 + κ / πR)

  e₁ = exp(-t)
  e₂ = exp(-(πR * α₂ + πY) * t)
  e₃ = exp(-(πY * α₁ + πR) * t)

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

  return Pmatrix(P₁,  P₅,  P₆,  P₅,
                 P₇,  P₂,  P₇,  P₈,
                 P₉,  P₁₀, P₃,  P₁₀,
                 P₁₁, P₁₂, P₁₁, P₄)
end
