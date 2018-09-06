abstract type JC69 <: NASM end


P(mod::JC69, t::Array{Float64}) = [P(mod, i) for i in t]


function distance(::Type{JC69}, seq1, seq2)
  match_matrix = _match_matrix(seq1, seq2)
  p = substitutions(match_matrix)
  l = sum(match_matrix)
  return Distance{JC69}(-(3.0/4.0) * log(1.0 - (4.0/3.0) * p), # mean
                        p * (1.0 - p) / ((1.0 - (4.0 * p / 3.0)^2) * l)) # variance
end


function distance(::Type{T}, seq1, seq2) where T <: JC69
  return distance(JC69, seq1, seq2)
end
