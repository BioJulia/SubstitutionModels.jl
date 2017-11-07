using Documenter, SubstitutionModels

makedocs(
    modules = [SubstitutionModels],
    format = :html,
    sitename = "SubstitutionModels.jl",
    pages = [
        "Home" => "index.md",
        "Manual" => [

        ],
        "Contributing" => "contributing.md"
    ],
    authors = "Justin Angevarre & Ben J. Ward, at the BioJulia organisation, and contributors."
)

deploydocs(
    deps = Deps.pip("mkdocs", "pygments", "mkdocs-material", "python-markdown-math"),
    repo = "github.com/BioJulia/SubstitutionModels.jl.git",
    julia = "0.6",
    osname = "linux",
)
