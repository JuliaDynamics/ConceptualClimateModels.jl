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
export construct_dynamical_system

include("constants.jl")
include("variables.jl")
include("statespace.jl")
include("default.jl")
include("processes_advanced.jl")

DEFAULT_DIFFEQ = DynamicalSystemsBase.DEFAULT_DIFFEQ

# TODO: Make in-place depend on state space dimension
function construct_dynamical_system(proc, default = DEFAULT_PROCESSES; diffeq = DEFAULT_DIFFEQ, inplace::Bool = true, kwargs...)
    sys = processes_to_mtkmodel(proc, default; kwargs...)
    ssys = structural_simplify(sys)
    # The usage of `nothing` for the initial state assumes all state variables
    # and all parameters have been defined with a default value
    if inplace
        prob = ODEProblem(ssys, nothing, (0.0, Inf))
    else
        prob = ODEProblem{false}(ssys, nothing, (0.0, Inf); u0_constructor = x->SVector(x...))
    end
    ds = CoupledODEs(prob, diffeq)
    return ds
end

end # module