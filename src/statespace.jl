# This function is only meaningful for dynamic variables!
"""
    plausible_limits(x::Num)

Return a tuple (min, max) of plausible limiting values for the variable `x`.
If `x` is defined with the `bounds` metadata, this is returned as the limits.
Else, if `x` has a default value, the limits are this value ± 20%.
Else, if there is no default value, an error is thrown.
"""
function plausible_limits(var)
    if ModelingToolkit.hasbounds(var)
        return ModelingToolkit.getbounds(var)
    elseif !isnothing(default_value(var))
        return (0.8default_value(var), 1.2default_value(var))
    else
        return throw(ArgumentError("""
        Variable $(var) does not have any defined limits.
        Define it with a `bounds` metadata or a default value.
        """))
    end
end

"""
    plausible_limits(ds::DynamicalSystem [, idxs])

Return a vector of limits (min, max) for each dynamic state variable in `ds`.
Optionally provide the `idxs` of the variables to use as a vector of `Symbol`s
for symbolic variables present in the referrenced MTK model of `ds`.

See also [`plausible_grid`](@ref), [`plausible_ic_sampler`](@ref).
"""
function plausible_limits(ds::DynamicalSystem, idxs = nothing)
    model = referrenced_sciml_model(ds)
    vars = ModelingToolkit.unknowns(model)
    if idxs isa Vector{Symbol}
        allvars = ModelingToolkit.unknowns(model)
        allnames = ModelingToolkit.getname.(allvars)
        j = [findfirst(isequal(i), allnames) for i in idxs]
        vars = allvars[j]
    end
    return plausible_limits.(vars)
end

"""
    plausible_ic_sampler(ds::DynamicalSystem [, seed])

Return a `sampler` that can be called as a 0-argument function `sampler()`, which
yields random initial conditions within the hyperrectangle defined by the
[`plausible_limits`](@ref) of `ds`.
The `sampler` is useful to give to e.g., `DynamicalSystems.basins_fractions`.
"""
function plausible_ic_sampler(ds::DynamicalSystem, seed = rand(1:10000))
    grid = plausible_grid(ds)
    sampler, = statespace_sampler(grid, seed)
    return sampler
end

"""
    plausible_grid(ds::DynamicalSystem, n = 101)

Return a `grid` that is a tuple of `range` objects that each spans the
[`plausible_limits`](@ref)(`ds`).
`n` can be an integer or a vector of integers (for each dimension).
The resulting `grid` is useful to give to `DynamicalSystems.AttractorsViaRecurrences`.
"""
function plausible_grid(ds::DynamicalSystem, n = 101)
    minmax = plausible_limits(ds)
    if n isa Integer
        ranges = [range(mm[1], mm[2]; length = n) for mm in minmax]
    elseif n isa AbstractVector
        ranges = [range(minmax[i][1], minmax[i][2]; length = n[i]) for i in eachindex(minmax)]
    end
    return Tuple(ranges)
end

export plausible_ic_sampler, plausible_limits, plausible_grid