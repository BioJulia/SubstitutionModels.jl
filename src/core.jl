abstract type SubstitutionModel end
const SM = SubstitutionModel

abstract type NucleicAcidSubstitutionModel <: SubstitutionModel end
const NASM = NucleicAcidSubstitutionModel

const Qmatrix = SMatrix{4, 4, Float64}
const Pmatrix = SMatrix{4, 4, Float64}
