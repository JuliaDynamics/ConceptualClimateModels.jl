export processes_to_coupledodes
export dynamical_system_summary
export all_equations
export named_current_parameters
export try_set_parameter!

DEFAULT_DIFFEQ = DynamicalSystemsBase.DEFAULT_DIFFEQ

"""
    processes_to_coupledodes(processes [, default]; kw...)

Convert a given `Vector` of processes to a `DynamicalSystem`, in particular `CoupledODEs`.
All processes represent symbolic equations managed by ModelingToolkit.jl.
`default` is a vector for default processes that "process-less" variables
introduced in `processes` will obtain.
Use [`processes_to_mtkmodel`](@ref) to obtain the MTK model before it is structurally
simplified and converted to a `DynamicalSystem`.
See also [`processes_to_mtkmodel`](@ref) for more details on what `processes` is,
or see the online [Tutorial](@ref).

## Keyword arguments

- `diffeq`: options passed to DifferentialEquations.jl ODE solving
  when constructing the `CoupledODEs`.
- `inplace`: whether the dynamical system will be in place or not.
  Defaults to `true` if the system dimension is ≤ 5.
- `split = false`: whether to split parameters as per ModelingToolkit.jl.
  Note the default is not ModelingToolkit's default, i.e., no splitting occurs.
  This accelerates parameter access, assuming all parameters are of the same type.
- `kw...`: all other keywords are propagated to `processes_to_mtkmodel`.
"""
function processes_to_coupledodes(proc, default = [];
        diffeq = DEFAULT_DIFFEQ, inplace::Bool = false, split::Bool = false, kwargs...
    )
    sys = processes_to_mtkmodel(proc, default; kwargs...)
    ssys = structural_simplify(sys; split)
    if isnothing(inplace)
        D = length(unknowns(ssys))
        inplace = D > 5
    end

    # The usage of `nothing` for the initial state assumes all state variables
    # and all parameters have been defined with a default value. This also means
    if inplace
        prob = ODEProblem(ssys, nothing, (0.0, Inf), nothing)
    else
        prob = ODEProblem{false}(ssys, nothing, (0.0, Inf), nothing; u0_constructor = x->SVector(x...))
    end
    ds = CoupledODEs(prob, diffeq)
    return ds
end

"""
    dynamical_system_summary(ds::DynamicalSystem)

Return a printable/writable string containing a summary of `ds`,
which outlines its current status and lists all symbolic
equations and parameters that constitute the system, if a referrence to a
ModelingToolkit.jl model exists in `ds`.
"""
function dynamical_system_summary(ebm::DynamicalSystem)
    if !DynamicalSystemsBase.has_referrenced_model(ebm)
        summary = sprint(show, MIME"text/plain"(), ebm)
        return summary
    end
    # Else, use symbolic stuff; first print the first 4 lines
    summary = keepfirstlines(sprint(show, MIME"text/plain"(), ebm), 4)*"\n"
    model = referrenced_sciml_model(ebm)
    summary *= "\nwith equations:\n"*skipfirstline(sprint(show, MIME"text/plain"(), all_equations(model)))
    # We can use
    # [`ModelingToolkit.dump_variable_metadata`](@ref), [`ModelingToolkit.dump_parameters`](@ref),
    # [`ModelingToolkit.dump_unknowns`](@ref), as well; but this works fine.
    summary *= "\n\nand parameter values:\n"*skipfirstline(sprint(show, MIME"text/plain"(), named_current_parameters(ebm)))
    return summary
end

keepfirstlines(str, n) = join(split(str, '\n')[1:n], '\n')
skipfirstline(str, limit = 2) = split(str, '\n'; limit)[2]

function ProcessBasedModelling.all_equations(ds::DynamicalSystem)
    model = referrenced_sciml_model(ds)
    return all_equations(model)
end

"""
    named_current_parameters(ds::DynamicalSystem)

Return a dictionary mapping parameters of `ds` (as `Symbol`s) to their values.
"""
function named_current_parameters(ds::DynamicalSystem)
    mtk = referrenced_sciml_model(ds)
    params_names = Symbol.(ModelingToolkit.parameters(mtk))
    params_values = current_parameter.(ds, params_names)
    params = Dict(params_names .=> params_values)
    return params
end

ProcessBasedModelling.has_symbolic_var(ds::DynamicalSystem, var) =
has_symbolic_var(referrenced_sciml_model(ds), var)

"""
    try_set_parameter!(ds::DynamicalSystem, symbol, value)

An extension of `set_parameter!` for as symbolic parameter,
that first checks if the symbol exists in the referenced MTK model,
and if not, it warns and returns `nothing`.
"""
function try_set_parameter!(ds::DynamicalSystem, symbol, value)
    if !has_symbolic_var(ds, symbol)
        @warn "Symbolic parameter $(p) does not exist in system"
        return nothing
    end
    set_parameter!(ds, symbol, value)
end
