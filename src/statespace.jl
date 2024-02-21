"""
    physically_plausible_limits(ds::DynamicalSystem)

Return a vector of limits (min, max) for each state variable in `ds`.
The states need to be named as the limits are deduced from the function
`physically_plausible_limits(varname::String)`.
"""
function physically_plausible_limits(ds::DynamicalSystem)
    model = referrenced_sciml_model(ds)
    vars = ModelingToolkit.states(model)
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