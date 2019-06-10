module SubstitutionModels

  using
    BioSymbols,
    StaticArrays,
    LinearAlgebra

  include("core.jl")
  include("indexing.jl")
  include("nucleic_acid/nucleic_acid.jl")

  for m in ["jc69"; "k80"; "f81"; "f84"; "hky85"; "tn93"; "gtr"]
    for f in ["abstract"; "absolute"; "relative"; "outer_constructors"]
      include("nucleic_acid/$(m)/$(f).jl")
    end
  end

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
