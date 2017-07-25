"Count the differences in two aligned DNA sequences, return the results in a
matrix in the form of:

  [[A→A, A→C, A→G, A→T]
   [C→A, C→C, C→G, C→T]
   [G→A, G→C, G→G, G→T]
   [T→A, T→C, T→G, T→T]]
"
function differences(seq1::BioSequence{DNAAlphabet{4}}, seq2::BioSequence{DNAAlphabet{4}})
  if length(seq1) != length(seq2)
    error("Sequences must be of the same length")
  end
  difference_matrix = MMatrix{4, 4, Int64}(0, 0, 0, 0,
                                           0, 0, 0, 0,
                                           0, 0, 0, 0,
                                           0, 0, 0, 0)
  @simd for i in 1:length(seq1)
    try
      difference_matrix[seq1[i], seq2[i]] += 1
    end
  end
  return difference_matrix
end

"Count the differences in two aligned RNA sequences, return the results in a
matrix in the form of:

  [[A→A, A→C, A→G, A→T]
   [C→A, C→C, C→G, C→T]
   [G→A, G→C, G→G, G→T]
   [T→A, T→C, T→G, T→T]]
"
function differences(seq1::BioSequence{RNAAlphabet{4}}, seq2::BioSequence{RNAAlphabet{4}})
  if length(seq1) != length(seq2)
    error("Sequences must be of the same length")
  end
  difference_matrix = MMatrix{4, 4, Int64}(0, 0, 0, 0,
                                           0, 0, 0, 0,
                                           0, 0, 0, 0,
                                           0, 0, 0, 0)
  @simd for i in 1:length(seq1)
    try
      difference_matrix[seq1[i], seq2[i]] += 1
    end
  end
  return difference_matrix
end

"Proportion of site showing differences. Referred to the observed distance or
p-distance"
function substitutions(d::MMatrix{4, 4, Int64})
  return 1.0 - (sum(diag(d)) / sum(d))
end

"Proportion of sites showing transitional differences. Referred to as the P
proportion"
function transitions(d::MMatrix{4, 4, Int64})
  return d[DNA_C, DNA_T] + d[DNA_T, DNA_C] +
         d[DNA_A, DNA_G] + d[DNA_G, DNA_A] /
         sum(d)
end

"Proportion of purine sites showing transitional differences. Referred to as the
 P1 proportion"
function purine_transitions(d::MMatrix{4, 4, Int64})
  return d[DNA_A, DNA_G] + d[DNA_G, DNA_A] /
         d[DNA_A, DNA_G] + d[DNA_G, DNA_A] +
         d[DNA_A, DNA_A] + d[DNA_G, DNA_G]
end

"Proportion of pyrimidine sites showing transitional differences. Referred to as
 the P2 proportion"
function pyrimidine_transitions(d::MMatrix{4, 4, Int64})
  return d[DNA_C, DNA_T] + d[DNA_T, DNA_C] /
         d[DNA_C, DNA_T] + d[DNA_T, DNA_C] +
         d[DNA_C, DNA_C] + d[DNA_T, DNA_T]
end

"Proportion of sites showing transversional differences. Referred to as the Q
proportion"
function transversions(d::MMatrix{4, 4, Int64})
  return d[DNA_A, DNA_C] + d[DNA_C, DNA_A] +
         d[DNA_G, DNA_T] + d[DNA_T, DNA_G] +
         d[DNA_A, DNA_T] + d[DNA_T, DNA_A] +
         d[DNA_C, DNA_G] + d[DNA_G, DNA_C] /
         sum(d)
end

substitutions(seq1, seq2) = substitutions(differences(seq1, seq2))
transitions(seq1, seq2) = transitions(differences(seq1, seq2))
purine_transitions(seq1, seq2) = purine_transitions(differences(seq1, seq2))
pyrimidine_transitions(seq1, seq2) = pyrimidine_transitions(differences(seq1, seq2))
transversions(seq1, seq2) = transversions(differences(seq1, seq2))
