cd(@__DIR__)

using EnergyBalanceModels

import Downloads
Downloads.download(
    "https://raw.githubusercontent.com/JuliaDynamics/doctheme/master/build_docs_with_style.jl",
    joinpath(@__DIR__, "build_docs_with_style.jl")
)
include("build_docs_with_style.jl")

pages =  [
    "Introduction" => "index.md",
    "Predefined processes" => "processes.md",
    "References" => "references.md",
]

using DocumenterCitations

bib = CitationBibliography(
    joinpath(@__DIR__, "references.bib");
    style=:authoryear
)

build_docs_with_style(pages, EnergyBalanceModels;
    authors = "George Datseris <datseris.george@gmail.com>",
    bib,
    # We need to remove the cross references because we don't list here
    # the whole `DynamicalSystem` API...
    warnonly = [:doctest, :missing_docs, :cross_references],
    repo = Remotes.GitHub("JuliaDynamics", "EnergyBalanceModels.jl")
)
