# TODO: Rename to sigmoid
# TODO: Logistic function is faster than tanh in Julia on windows 10
export TanhProcess

"""
    TanhProcess(variable, driver, left, right, scale, reference) <: Process
    TanhProcess(variable, driver; left, right, scale, reference) <: Process
    TanhProcess(variable, driver; left, right, scale, start) <: Process

A common process for when a `variable` has a tanh-dependence on a `driver` variable.
The rest of the input arguments should be real numbers or `@parameter` named parameters.

The process creates the expression:
```
variable ~ left + (right - left)*(1 + tanh(2(driver - reference)/scale))*0.5
```
i.e., a tanh formula that goes from value `left` to value `right` as a function of `driver`
over a range of `scale` being centered at `reference`.
Alternatively, you can provide `start` instead of reference,
which creates the reference value as `reference = start - scale/2`.

If the values given to the parameters of the expression are real numbers, they become
named parameters prefixed with the name of `variable`, then the name of the `driver`,
and then `_tanh_left`, `_tanh_right`, `_tanh_rate` and `_tanh_ref` respectively.
Use `LiteralParameter` for parameters you do not wish to rename.
"""
struct TanhProcess <: Process
    variable
    driver_variable
    left
    right
    scale
    reference
end

function TanhProcess(variable, driver; left=0, right=1, scale=1, reference=nothing, start=nothing)
    if !isnothing(reference) && !isnothing(start)
        error("Only one of `reference, start` can be given!")
    end
    if isnothing(reference)
        reference = start - scale/2
    end
    return TanhProcess(variable, driver, left, right, scale, reference)
end

function ProcessBasedModelling.rhs(p::TanhProcess)
    y = p.variable
    x = p.driver_variable
    # Create the names for everything
    (; left, right, scale, reference) = p
    values = left, right, scale, reference
    xname = string(ModelingToolkit.getname(x))
    suffixes = (xname*"_") .* ["tanh_left", "tanh_right", "tanh_scale", "tanh_ref"]
    mtk_vars = map((val, s) -> new_derived_named_parameter(y, val, s; prefix = false), values, suffixes)
    # pass them to the regular tanh Julia function
    return tanh_expression(x, mtk_vars...)
end
function tanh_expression(T, left, right, scale, reference)
    return left + (right - left)*(1 + tanh(2(T - reference)/(scale)))*0.5
end
