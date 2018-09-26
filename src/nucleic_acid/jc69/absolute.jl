mutable struct JC69abs <: JC69
    λ::Float64
    function JC69abs(λ::Float64)
        if λ <= 0.
            error("JC69 parameter λ must be positive")
        end
        new(λ)
    end
end

"""
```julia-repl
julia> model = JC69(1.0);

julia> setrate!(model, [0.5])
Jukes and Cantor 1969 model (absolute rate form)
λ = 0.5
```
"""
@inline function setrate!(mod::JC69abs, rate::Vector{Float64})
    #check length of vector
    if length(rate) != 1
      @error "JC69 rate must be a vector of length 1"
    else 
      #separate vector into pieces
      λ = rate[1]
      #check correctness of rates (from above)
      if λ <= 0.
        @error "JC69 parameter λ must be positive"
      end
      #add rate
      mod.λ = λ
    end
    return mod
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
  P₂ = 0.25 + 0.25 * ω
  return Pmatrix(P₁, P₂, P₂, P₂,
                 P₂, P₁, P₂, P₂,
                 P₂, P₂, P₁, P₂,
                 P₂, P₂, P₂, P₁)
end
