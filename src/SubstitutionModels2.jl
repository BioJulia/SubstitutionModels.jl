__precompile__()

module SubstitutionModels2

  using
    BioSymbols,
    StaticArrays

  import
    Base.getindex

  include("core.jl")
  include("indexing.jl")
  include("nucleic_acid/jc69.jl")
  include("nucleic_acid/nucleic_acid.jl")

  export
    SubstitutionModel,
    SM,
    NucleicAcidSubstitutionModel,
    NASM,
    # AminoAcidSubstitutionModel,
    JC69,
    JC69abs,
    JC69rel,
    Q,
    P,
    π,
    μ

end # module
