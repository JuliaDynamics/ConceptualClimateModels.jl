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

include("processes_advanced.jl")
include("statespace.jl")
include("dynamicalsystems.jl")

include("GlobalMeanEBM/GlobalMeanEBM.jl")
include("CloudToppedMixedLayerModel/CloudToppedMixedLayerModel.jl")

end # module