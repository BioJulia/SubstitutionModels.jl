__precompile__()

module SubstitutionModels

  using
    BioSymbols,
    BioSequences,
    StaticArrays

  import
    Base.getindex,
    Base.setindex!,
    Base.show

  include("core.jl")
  include("indexing.jl")
  include("differences.jl")
  include("nucleic_acid/jc69.jl")
  include("nucleic_acid/k80.jl")
  include("nucleic_acid/f81.jl")
  include("nucleic_acid/f84.jl")
  include("nucleic_acid/hky85.jl")
  include("nucleic_acid/tn93.jl")
  include("nucleic_acid/gtr.jl")
  include("nucleic_acid/nucleic_acid.jl")

  export
    SubstitutionModel, SM,
    NucleicAcidSubstitutionModel, NASM,
    differences,
    JC69, JC69abs, JC69rel,
    K80, K80abs, K80rel,
    F81, F81abs, F81rel,
    F84, F84abs, F84rel,
    HKY85, HKY85abs, HKY85rel,
    TN93, TN93abs, TN93rel,
    GTR, GTRabs, GTRrel,
    Q, P, P_generic,
    _π

end # module
