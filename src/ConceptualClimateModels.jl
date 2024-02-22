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
export processes_to_coupledodes

include("constants.jl")
include("variables.jl")
include("statespace.jl")
include("default.jl")
include("processes_advanced.jl")

DEFAULT_DIFFEQ = DynamicalSystemsBase.DEFAULT_DIFFEQ

"""
    processes_to_coupledodes(processes, default = DEFAULT_CCM_PROCESSES; kw...)

Convert a given `Vector` of processes to a `DynamicalSystem`, in particular `CoupledODEs`.
All processes represent symbolic equations managed by ModelingToolkit.jl.
`default` is a vector for default processes that "process-less" variables
introduced in `processes` will obtain.
Use [`process_to_mtkmodel`](@ref) to obtain the MTK model before it is structurally
simplified and converted to a `DynamicalSystem`.
See also [`process_to_mtkmodel`](@ref) for more details on what `processes` is,
or see the online [Tutorial](@ref).

## Keyword arguments

- `diffeq`: options passed to DifferentialEquations.jl ODE solving
  when constructing the `CoupledODEs`.
- `inplace`: whether the dynamical system will be in place or not.
  Defaults to `true` if the system dimension is ≤ 5.
- `split = false`: whether to split parameters as per ModelingToolkit.jl.
  Note the default is not ModelingToolkit's default, i.e., no splitting occurs.
  This accelerates parameter access, assuming all parameters are of the same type.
- `kw...`: all other keywords are propagated to `process_to_mtkmodel`.
"""
function processes_to_coupledodes(proc, default = DEFAULT_CCM_PROCESSES;
        diffeq = DEFAULT_DIFFEQ, inplace = nothing, split::Bool = false, kwargs...
    )
    sys = processes_to_mtkmodel(proc, default; kwargs...)
    ssys = structural_simplify(sys; split)
    if isnothing(inplace)
        D = length(unknowns(ssys))
        inplace = D < 6
    end

    # The usage of `nothing` for the initial state assumes all state variables
    # and all parameters have been defined with a default value. This also means
    # that we can use `default_value` to
    if inplace
        prob = ODEProblem(ssys, nothing, (0.0, Inf))
    else
        prob = ODEProblem{false}(ssys, nothing, (0.0, Inf); u0_constructor = x->SVector(x...))
    end
    ds = CoupledODEs(prob, diffeq)
    return ds
end

end # module