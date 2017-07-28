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


const K80(α, β) = K80abs(α, β)


@inline function Q(mod::K80abs)
  α = mod.α; β = mod.β

  Q₁ = 0.25 * α
  Q₂ = 0.25 * β
  Q₃ = -(2 * Q₂ + Q₁)

  return Qmatrix(Q₃, Q₂, Q₁, Q₂,
                 Q₂, Q₃, Q₂, Q₁,
                 Q₁, Q₂, Q₃, Q₂,
                 Q₂, Q₁, Q₂, Q₃)
end


@inline function P(mod::K80abs, t::Float64)
  if t < 0.0
    error("t must be positive")
  end
  α = mod.α; β = mod.β

  e₁ = 0.25 * exp(-4 * β * t)
  e₂ = 0.5 * exp(-2 * (α + β) * t)

  P₁ = 0.25 + e₁ + e₂
  P₂ = 0.25 + e₁ - e₂
  P₃ = 0.25 - e₁

  return Pmatrix(P₁, P₃, P₂, P₃,
                 P₃, P₁, P₃, P₂,
                 P₂, P₃, P₁, P₃,
                 P₃, P₂, P₃, P₁)
end
