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
  include("nucleic_acid/nucleic_acid.jl")

  export
    SubstitutionModel,
    SM,
    NucleicAcidSubstitutionModel,
    NASM,
    JC69,
    JC69abs,
    JC69rel,
    K80,
    K80abs,
    K80rel,
    F81,
    F81abs,
    F81rel,
    Q, P,
    _π,
    _πR, _πY,
    _πA, _πC, _πG, _πT,
    _μ,
    _α, _β, _γ, _δ, _ϵ, _η

end # module
