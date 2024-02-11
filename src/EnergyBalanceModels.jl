module EnergyBalanceModels

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end EnergyBalanceModels


# Rules: all variables have a default value.
# All processes are represented via `Equation`s. Everything returns `Equation`s.
# Specifications are given as Spec files to the main equation.
# TODO:
# The spec files specify what processes to use for each fundamental component
# of the energy balance equations. Processes that are not given
# specs are not included (or have a dummy default spec that makes them be
# constants; e.g., albedo of clouds being 0.1 or so.)
# Multiple dispatch is then used to create the equations for each processes
# and also to decide whether to include any additional processes.
# We can find the spec files by their subtyping relation

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