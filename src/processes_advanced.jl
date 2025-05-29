# TODO: Logistic function is faster than tanh in Julia on windows 10
# Is it worth re-writing the Sigmoid process on logistic function?
export SigmoidProcess, ClampedLinearProcess

"""
    SigmoidProcess <: Process
    SigmoidProcess(variable, driver; left, right, scale, reference)

A common process for when a `variable` has a sigmoidal dependence on a `driver` variable.
The rest of the input arguments should be real numbers or `@parameter` named parameters.

The process creates a sigmoidal relationship based on the `tanh` function:
```
variable ~ left + (right - left)*(1 + tanh(2(driver - reference)/scale))*0.5
```
i.e., the variable goes from value `left` to value `right` as `driver` increases
over a range of `scale` (centered at `reference`).
Instead of `reference` you may provide `start` or `finish` keywords,
which make `reference = start + scale/2` or `reference = finish - scale/2` respectively.

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

function SigmoidProcess(variable, driver; left=0, right=1, scale=1,
    reference=nothing, start=nothing, finish = nothing)
    if count(!isnothing, (reference, start, finish)) > 1
        error("Exactly one of `reference, start, finish` must be given!")
    end
    if !isnothing(start)
        reference = start + scale/2
    elseif !isnothing(finish)
        reference = finish - scale/2
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

"""
    ClampedLinearProcess <: Process
    ClampedLinearProcess(variable, driver; left, right, left_driver, right_driver)

A common process for when a `variable` has a clamped linear dependence on a `driver` variable.
It is similar with the `SigmoidalProcess` but not smooth.

The `variable` increases linearly as a function of `driver`.
It increases from value `left` at driver's value `left_driver` to value
`right` at driver's value `right_driver`. The `variable` stays
at its left/right value respectively when the driver exceeds the given bounds.

The input keywords can all be real numbers or `@parameter` named parameters.
"""
struct ClampedLinearProcess <: Process
    variable
    driver_variable
    left
    right
    left_driver
    right_driver
end

function ClampedLinearProcess(variable, driver_variable; left = 0, right = 1, left_driver = 1, right_driver = 1)
    return ClampedLinearProcess(variable, driver_variable, left, right, left_driver, right_driver)
end

function ProcessBasedModelling.rhs(p::ClampedLinearProcess)
    x = p.driver_variable
    # Create the names for everything
    (; left, right, left_driver, right_driver) = p
    expr = left + (x - left_driver)*(right - left)/(right_driver - left_driver)
    expr = clamp(expr, left, right)
    return expr
end