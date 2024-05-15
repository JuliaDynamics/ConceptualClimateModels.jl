cd(@__DIR__)

using ConceptualClimateModels

import Downloads
Downloads.download(
    "https://raw.githubusercontent.com/JuliaDynamics/doctheme/master/build_docs_with_style.jl",
    joinpath(@__DIR__, "build_docs_with_style.jl")
)
include("build_docs_with_style.jl")

pages =  [
    "Introduction" => "index.md",
    "Tutorial" => "tutorial.md",
    "Submodules" => [
        "Global mean EBM" => "globalmeanebm.md",
    ],
    "Examples" => "examples.md",
    "References" => "references.md",
]

using DocumenterCitations

bib = CitationBibliography(
    joinpath(@__DIR__, "references.bib");
    style=:authoryear
)

build_docs_with_style(pages, ConceptualClimateModels, ProcessBasedModelling;
    authors = "George Datseris <datseris.george@gmail.com>",
    bib, warnonly = [:doctest, :missing_docs, :cross_references],
    repo = Remotes.GitHub("JuliaDynamics", "ConceptualClimateModels.jl")
)
