abstract type SubstitutionModel end


const SM = SubstitutionModel


"""
`NucleicAcidSubstitutionModel` is an abstract type that contains all models
describing a substitution process impacting biological sequences of `DNA` or
`RNA` with continous time Markov models.
"""
abstract type NucleicAcidSubstitutionModel <: SubstitutionModel end


const NASM = NucleicAcidSubstitutionModel


const Qmatrix = SMatrix{4, 4, Float64}


const Pmatrix = SMatrix{4, 4, Float64}
