module ConceptualClimateModels

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end ConceptualClimateModels

using Reexport
@reexport using ProcessBasedModelling # exports `t` as time
@reexport using DynamicalSystemsBase

include("constants.jl")
include("variables.jl")
include("statespace.jl")
include("default.jl")
include("processes_advanced.jl")

end # module