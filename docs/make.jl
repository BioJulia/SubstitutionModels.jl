using Documenter, SubstitutionModels

makedocs()

deploydocs(
    deps = Deps.pip("mkdocs", "pygments", "mkdocs-material", "python-markdown-math"),
    repo = "github.com/BioJulia/SubstitutionModels.jl.git",
    julia = "0.6",
    osname = "linux",
)
