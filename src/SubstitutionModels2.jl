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
  include("nucleic_acid/k80.jl")
  include("nucleic_acid/f81.jl")
  include("nucleic_acid/f84.jl")
  include("nucleic_acid/hky85.jl")
  include("nucleic_acid/tn93.jl")
  include("nucleic_acid/nucleic_acid.jl")

  export
    SubstitutionModel, SM,
    NucleicAcidSubstitutionModel, NASM,
    JC69, JC69abs, JC69rel,
    K80, K80abs, K80rel,
    F81, F81abs, F81rel,
    F84, F84abs, F84rel,
    HKY85, HKY85abs, HKY85rel,
    TN93, TN93abs, TN93rel,
    Q, P,
    _π

end # module
