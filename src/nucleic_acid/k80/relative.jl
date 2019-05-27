struct K80rel <: K80
  κ::Float64
  function K80rel(κ::Float64)
    if κ <= 0.
      error("K80 parameter κ must be positive")
    end
    new(κ)
  end
end


function show(io::IO, object::K80rel)
  print(io, "\r\e[0m\e[1mK\e[0mimura 19\e[1m80\e[0m model (relative rate form)
κ = $(object.κ)")
end


const K80(κ) = K80rel(κ)


@inline function Q(mod::K80rel)
  κ = mod.κ

  Q₁ = 0.25 * κ
  Q₂ = 0.25
  Q₃ = -(2 * Q₂ + Q₁)

  return Qmatrix(Q₃, Q₂, Q₁, Q₂,
                 Q₂, Q₃, Q₂, Q₁,
                 Q₁, Q₂, Q₃, Q₂,
                 Q₂, Q₁, Q₂, Q₃)
end


@inline function P(mod::K80rel, t::Float64)
  if t < 0.0
    error("t must be positive")
  end
  κ = mod.κ

  e₁ = 0.25 * exp(-t)
  e₂ = 0.5 * exp(-0.5 * t * (κ + 1))

  P₁ = 0.25 + e₁ + e₂
  P₂ = 0.25 + e₁ - e₂
  P₃ = 0.25 - e₁

  return Pmatrix(P₁, P₃, P₂, P₃,
                 P₃, P₁, P₃, P₂,
                 P₂, P₃, P₁, P₃,
                 P₃, P₂, P₃, P₁)
end
