struct JC69abs <: JC69
    λ::Float64
    function JC69abs(λ::Float64)
        if λ <= 0.
            error("JC69 parameter λ must be positive")
        end
        new(λ)
    end
end


function show(io::IO, object::JC69abs)
  print(io, "\r\e[0m\e[1mJ\e[0mukes and \e[1mC\e[0mantor 19\e[1m69\e[0m model (absolute rate form)
λ = $(object.λ)")
end


JC69(λ) = JC69abs(λ)


function Q(mod::JC69abs)
  λ = mod.λ
  Q₁ =  0.25 * λ
  Q₂ = -(Q₁ * 3)
  return Qmatrix(Q₂, Q₁, Q₁, Q₁,
                 Q₁, Q₂, Q₁, Q₁,
                 Q₁, Q₁, Q₂, Q₁,
                 Q₁, Q₁, Q₁, Q₂)
end


function P(mod::JC69abs, t::Float64)
  if t < 0.0
    error("t must be positive")
  end
  λ = mod.λ
  ω = exp(-t * λ)
  P₁ = 0.25 + 0.75 * ω
  P₂ = 0.25 - 0.25 * ω
  return Pmatrix(P₁, P₂, P₂, P₂,
                 P₂, P₁, P₂, P₂,
                 P₂, P₂, P₁, P₂,
                 P₂, P₂, P₂, P₁)
end
