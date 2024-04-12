# This function is only meaningful for dynamic variables!
"""
    physically_plausible_limits(x)

Return a tuple (min, max) of plausible limiting values for the variable `x`.
If the variable does not have defined `bounds` metadata, then the default value ± 20% is used.
If there is no default value, a heuristic is tried, and an error is thrown if it fails.
"""
function physically_plausible_limits(var)
    if ModelingToolkit.hasbounds(var)
        return getbounds(var)
    elseif !isnothing(default_value(var))
        return (0.8default_value(var), 1.2default_value(var))
    else
        return CCMV.physically_plausible_limits(string(ModelingToolkit.getname(var)))
    end
end

"""
    physically_plausible_limits(ds::DynamicalSystem [, idxs])

Return a vector of limits (min, max) for each dynamic state variable in `ds`.
Optionally provide the `idxs` of the variables to use as a vector of `Symbol`s
or a vector of symbolic variables present in the referrenced MTK model of `ds`.
"""
function physically_plausible_limits(ds::DynamicalSystem, idxs = nothing)
    model = referrenced_sciml_model(ds)
    vars = ModelingToolkit.unknowns(model)
    if idxs isa Nothing
        vars = ModelingToolkit.unknowns(model)
    elseif idxs isa Vector{Symbol}
        vars = ModelingToolkit.unknowns(model)
        vars = filter!(v -> ModelingToolkit.getname(v) ∈ idxs, vars)
    else
        vars = idxs
    end
    minmax = physically_plausible_limits.(vars)
    return minmax
end

"""
    physically_plausible_ic_sampler(ds::DynamicalSystem)

Return a `sampler` that can be called as a 0-argument function `sampler()`, which
yields random initial conditions within the hyperrectangle defined by the
[`physically_plausible_limits`](@ref) of `ds`.
The `sampler` is useful to give to e.g., `DynamicalSystems.basins_fractions`.
"""
function physically_plausible_ic_sampler(ds::DynamicalSystem)
    minmax = physically_plausible_limits(ds)
    mini = getindex.(minmax, 1)
    maxi = getindex.(minmax, 2)
    sampler, = statespace_sampler(HRectangle(mini, maxi))
    return sampler
end

"""
    physically_plausible_grid(ds::DynamicalSystem, n = 101)

Return a `grid` that is a tuple of `range` objects that each spans the
[`physically_plausible_limits`](@ref) of `ds`.
`n` can be an integer or a vector of integers (for each dimension).
The resulting `grid` is useful to give to `DynamicalSystems.AttractorsViaRecurrences`.
"""
function physically_plausible_grid(ds::DynamicalSystem, n = 101)
    minmax = physically_plausible_limits(ds)
    if n isa Integer
        ranges = [range(mm[1], mm[2]; length = n) for mm in minmax]
    elseif n isa AbstractVector
        ranges = [range(minmax[i][1], minmax[i][2]; length = n[i]) for i in eachindex(minmax)]
    end
    return Tuple(ranges)
end

export physically_plausible_ic_sampler,
    physically_plausible_limits,
    physically_plausible_grid