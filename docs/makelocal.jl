using Documenter, SubstitutionModels

makedocs(
format = :html,
sitename = "SubstitutionModels.jl",
pages = ["index.md"]
)

deploydocs(
    deps = Deps.pip("mkdocs", "pygments", "mkdocs-material", "python-markdown-math"),
    repo = "github.com/BioJulia/SubstitutionModels.jl.git",
    julia = "0.6",
    osname = "linux",
)
