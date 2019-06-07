using Documenter, SubstitutionModels

makedocs(
    modules = [SubstitutionModels],
    format = Documenter.HTML(),
    sitename = "SubstitutionModels.jl",
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Substitution models" => "man/models.md",
            "Q and P matrices" => "man/pqmatrices.md",
            "Provided models & custom models" => "man/modeltypes.md"
        ],
        "Contributing" => "contributing.md"
    ],
    authors = "Justin Angevarre & Ben J. Ward, at the BioJulia organisation, and contributors."
)

deploydocs(
    repo = "github.com/BioJulia/SubstitutionModels.jl.git"
)
