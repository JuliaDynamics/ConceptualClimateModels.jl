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

# include all processes first
for (root, dirs, files) in walkdir(joinpath(@__DIR__, "processes"))
    for file in files
        include(joinpath(root, file))
    end
end

include("default.jl")
export DEFAULT_CCM_PROCESSES

include("processes_advanced.jl")
include("dynamicalsystems.jl")

export DEFAULT_PROCESSES

end # module