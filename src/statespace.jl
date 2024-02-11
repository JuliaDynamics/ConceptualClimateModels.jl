"""
    physically_plausible_limits(ds::DynamicalSystem)

Return a vector of limits (min, max) for each state variable in `ds`.
The states need to be named as the limits are deduced from the function
`physically_plausible_limits(varname::String)`.
"""
function physically_plausible_limits(ds::DynamicalSystem)
    sys = referrenced_sciml_model(ds)
    vars = ModelingToolkit.states(sys)
    varstrings = @. string(ModelingToolkit.getname(vars))
    minmax = physically_plausible_limits.(varstrings)
    return minmax
end


function physically_plausible_ic_sampler(ds::DynamicalSystem)
    minmax = physically_plausible_limits(ds)
    mini = getindex.(minmax, 1)
    maxi = getindex.(minmax, 2)
    sampler, = statespace_sampler(HRectangle(mini, maxi))
    return sampler
end

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