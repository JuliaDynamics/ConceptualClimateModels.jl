export processes_to_coupledodes
export dynamical_system_summary

DEFAULT_DIFFEQ = DynamicalSystemsBase.DEFAULT_DIFFEQ

"""
    processes_to_coupledodes(processes, default = DEFAULT_PROCESSES; kw...)

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
function processes_to_coupledodes(proc, default = DEFAULT_PROCESSES;
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

"""
    dynamical_system_summary(ds::DynamicalSystem)

Return a printable string containing a summary of `ds`,
which outlines its current status and lists all symbolic
equations and parameters that constitute the system, if a referrence to a
ModelingToolkit.jl exists in `ds`.
"""
function dynamical_system_summary(ebm::DynamicalSystem)
    if !DynamicalSystemsBase.has_referrenced_model(ebm)
        summary = sprint(show, MIME"text/plain"(), ebm)
        return summary
    end
    # Else, use symbolic stuff
    summary = keepfirstlines(sprint(show, MIME"text/plain"(), ebm), 4)*"\n"
    model = referrenced_sciml_model(ebm)
    summary *= "\nwith dynamic equations:\n"*skipfirstline(sprint(show, MIME"text/plain"(), equations(model)))
    summary *= "\n\nwith observed variables:\n"*skipfirstline(sprint(show, MIME"text/plain"(), observed(model)))
    summary *= "\n\nand parameter values:\n"*skipfirstline(sprint(show, MIME"text/plain"(), named_current_parameters(ebm)))
    return summary
end

keepfirstlines(str, n) = join(split(str, '\n')[1:n], '\n')
skipfirstline(str, limit = 2) = split(str, '\n'; limit)[2]

function named_current_parameters(ebm)
    cp = current_parameters(ebm)
    pnames = parameters(referrenced_sciml_model(ebm))
    return [pn ~ p for (p, pn) in zip(cp, pnames)]
end
