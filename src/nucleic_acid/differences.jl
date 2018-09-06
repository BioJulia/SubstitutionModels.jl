function _match_matrix(seq1::T, seq2::T) where T <: Union{BioSequence{DNAAlphabet{4}}, BioSequence{RNAAlphabet{4}}}
  if length(seq1) != length(seq2)
    @error "Sequences must be of the same length"
  end
  ambiguities = 0
  counts = MMatrix{4, 4, Int64}(0, 0, 0, 0,
                                0, 0, 0, 0,
                                0, 0, 0, 0,
                                0, 0, 0, 0)
  @simd for i in 1:length(seq1)
    try
      counts[seq1[i], seq2[i]] += 1
    catch
      ambiguities += 1
    end
  end
  if ambiguities > 0
    @info "$ambiguities ambiguious comparisions excluded"
  end
  return counts
end


"""
F matrix - A matrix of `NucleicAcid` proportions of aligned `BioSequences`
"""
function F_hat(d::T) where T<:Union{MMatrix{4, 4, Int64}, SMatrix{4, 4, Int64}}
  return d ./ sum(d)
end


"""
A symmetricized F matrix
"""
function F_hat_symmetricized(d)
  f = F_hat(d)
  return (f + f')/2.0
end


"""Proportion of sites showing differences. Referred to the observed distance or p-distance
"""
function substitutions(d)
  return 1.0 - sum(diag(F_hat(d)))
end


"""
Proportion of sites showing transitional differences. Referred to as the P proportion
"""
function transitions(d)
  f = F_hat(d)
  return f[DNA_C, DNA_T] + f[DNA_T, DNA_C] + f[DNA_A, DNA_G] + f[DNA_G, DNA_A]
end


"""
Proportion of purine sites showing transitional differences. Referred to as the P1 proportion
"""
function purine_transitions(d::T) where T<:Union{MMatrix{4, 4, Int64}, SMatrix{4, 4, Int64}}
  return d[DNA_A, DNA_G] + d[DNA_G, DNA_A] /
         d[DNA_A, DNA_G] + d[DNA_G, DNA_A] +
         d[DNA_A, DNA_A] + d[DNA_G, DNA_G]
end


"""
Proportion of pyrimidine sites showing transitional differences. Referred to as the P2 proportion
"""
function pyrimidine_transitions(d::T) where T<:Union{MMatrix{4, 4, Int64}, SMatrix{4, 4, Int64}}
  return d[DNA_C, DNA_T] + d[DNA_T, DNA_C] /
         d[DNA_C, DNA_T] + d[DNA_T, DNA_C] +
         d[DNA_C, DNA_C] + d[DNA_T, DNA_T]
end


"""
Proportion of sites showing transversional differences. Referred to as the Q proportion
"""
function transversions(d)
  f = F_hat(d)
  return f[DNA_A, DNA_C] + f[DNA_C, DNA_A] +
         f[DNA_G, DNA_T] + f[DNA_T, DNA_G] +
         f[DNA_A, DNA_T] + f[DNA_T, DNA_A] +
         f[DNA_C, DNA_G] + f[DNA_G, DNA_C]
end


F_hat(seq1, seq2) = F_hat(_match_matrix(seq1, seq2))
F_hat_symmetricized(seq1, seq2) = F_hat_symmetricized(_match_matrix(seq1, seq2))
substitutions(seq1, seq2) = substitutions(_match_matrix(seq1, seq2))
transitions(seq1, seq2) = transitions(_match_matrix(seq1, seq2))
purine_transitions(seq1, seq2) = purine_transitions(_match_matrix(seq1, seq2))
pyrimidine_transitions(seq1, seq2) = pyrimidine_transitions(_match_matrix(seq1, seq2))
transversions(seq1, seq2) = transversions(_match_matrix(seq1, seq2))
