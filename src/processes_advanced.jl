# TODO: Logistic function is faster than tanh in Julia on windows 10
export SigmoidProcess

"""
    SigmoidProcess <: Process
    SigmoidProcess(variable, driver; left, right, scale, reference)
    SigmoidProcess(variable, driver; left, right, scale, start)

A common process for when a `variable` has a sigmoidal dependence on a `driver` variable.
The rest of the input arguments should be real numbers or `@parameter` named parameters.

The process creates a sigmoidal relationship based on the `tanh` function:
```
variable ~ left + (right - left)*(1 + tanh(2(driver - reference)/scale))*0.5
```
i.e., the variable goes from value `left` to value `right` as `driver` increases
over a range of `scale` (centered at `reference`).
Alternatively, you can provide `start` instead of `reference`,
which creates the reference value as `reference = start - scale/2`.

If the values given to the parameters of the expression are real numbers, they become
named parameters prefixed with the name of `variable`, then the name of the `driver`,
and then `_sigmoid_left`, `_sigmoid_right`, `_sigmoid_rate` and `_sigmoid_ref` respectively.
Use `LiteralParameter` for parameters you do not wish to rename.
"""
struct SigmoidProcess <: Process
    variable
    driver_variable
    left
    right
    scale
    reference
end

function SigmoidProcess(variable, driver; left=0, right=1, scale=1, reference=nothing, start=nothing)
    if !isnothing(reference) && !isnothing(start)
        error("Only one of `reference, start` can be given!")
    end
    if isnothing(reference)
        reference = start - scale/2
    end
    return SigmoidProcess(variable, driver, left, right, scale, reference)
end

function ProcessBasedModelling.rhs(p::SigmoidProcess)
    y = p.variable
    x = p.driver_variable
    # Create the names for everything
    (; left, right, scale, reference) = p
    values = left, right, scale, reference
    xname = string(ModelingToolkit.getname(x))
    suffixes = (xname*"_") .* ["sigmoid_left", "sigmoid_right", "sigmoid_scale", "sigmoid_ref"]
    mtk_vars = map((val, s) -> new_derived_named_parameter(y, val, s; prefix = false), values, suffixes)
    # pass them to the regular tanh Julia function
    return sigmoid_expression(x, mtk_vars...)
end
function sigmoid_expression(T, left, right, scale, reference)
    return left + (right - left)*(1 + tanh(2(T - reference)/(scale)))*0.5
end
