module SubstitutionModels

  using
    BioSymbols,
    StaticArrays,
    LinearAlgebra

  import
    Base.getindex,
    Base.show,
    Base.setindex!,
    Base.checkbounds

  include("core.jl")
  include("indexing.jl")
  include("nucleic_acid/nucleic_acid.jl")
  include("nucleic_acid/jc69/abstract.jl")
  include("nucleic_acid/jc69/absolute.jl")
  include("nucleic_acid/jc69/relative.jl")
  include("nucleic_acid/k80/abstract.jl")
  include("nucleic_acid/k80/absolute.jl")
  include("nucleic_acid/k80/relative.jl")
  include("nucleic_acid/f81/abstract.jl")
  include("nucleic_acid/f81/absolute.jl")
  include("nucleic_acid/f81/relative.jl")
  include("nucleic_acid/f84/abstract.jl")
  include("nucleic_acid/f84/absolute.jl")
  include("nucleic_acid/f84/relative.jl")
  include("nucleic_acid/hky85/abstract.jl")
  include("nucleic_acid/hky85/absolute.jl")
  include("nucleic_acid/hky85/relative.jl")
  include("nucleic_acid/tn93/abstract.jl")
  include("nucleic_acid/tn93/absolute.jl")
  include("nucleic_acid/tn93/relative.jl")
  include("nucleic_acid/gtr/abstract.jl")
  include("nucleic_acid/gtr/absolute.jl")
  include("nucleic_acid/gtr/relative.jl")

  export
    SubstitutionModel, SM,
    NucleicAcidSubstitutionModel, NASM,
    JC69, JC69abs, JC69rel,
    K80, K80abs, K80rel,
    F81, F81abs, F81rel,
    F84, F84abs, F84rel,
    HKY85, HKY85abs, HKY85rel,
    TN93, TN93abs, TN93rel,
    GTR, GTRabs, GTRrel,
    Q, P

end # module
