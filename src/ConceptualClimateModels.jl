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
using ProcessBasedModelling: rhs, lhs_variable, lhs

import NaNMath # so we can use them for sqrt, log, etc. in MTK

include("constants.jl")

include("variables.jl")
using .CCMV

include("statespace.jl")

export DEFAULT_CCM_PROCESSES

include("processes_advanced.jl")
include("dynamicalsystems.jl")

export DEFAULT_PROCESSES

end # module