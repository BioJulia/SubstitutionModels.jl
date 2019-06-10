struct JC69rel <: JC69 end


function Base.show(io::IO, object::JC69rel)
  print(io, "\r\e[0m\e[1mJ\e[0mukes and \e[1mC\e[0mantor 19\e[1m69\e[0m model (relative rate form)")
end


@inline function Q(mod::JC69rel)
  Q₁ =  0.25
  Q₂ = -(Q₁ * 3)
  return SMatrix{4, 4, Float64}(Q₂, Q₁, Q₁, Q₁,
                                Q₁, Q₂, Q₁, Q₁,
                                Q₁, Q₁, Q₂, Q₁,
                                Q₁, Q₁, Q₁, Q₂)
end


@inline function P(mod::JC69rel, t::Float64)
  if t < 0.0
    error("t must be positive")
  end
  e₁ = exp(-t)
  P₁ = 0.25 + 0.75 * e₁
  P₂ = 0.25 - 0.25 * e₁
  return Pmatrix(P₁, P₂, P₂, P₂,
                 P₂, P₁, P₂, P₂,
                 P₂, P₂, P₁, P₂,
                 P₂, P₂, P₂, P₁)
end
