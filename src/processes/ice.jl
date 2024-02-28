export IceAlbedoFeedback

# All of this could be made a simple function that returns the
# ExpRelaxation(Tanh(..)) but having it this way makes for better printing
# and more appropriately-named parameters.
"""
    IceAlbedoFeedback(; T, α_ice,
        max = 0.45, min = 0.1, Tscale = 10, Tfreeze = 275.15, τ = 0
    )

Create an equation that assigns ice albedo `α_ice` to a hyperbolic tangent of temperature `T`.
This represents an approximately linear decrease with `T`, as ice melts over part of the
earth, while it is constant for all `T` for which the earth would be either entirely ice
covered (`T < Tfreeze - scale`) or ice free (`T > Tfreeze`).

In essence this is a [`TanhProcess`](@ref) with the given keywords as parameters
with reference temperature `Tref = Tfreeze - scale/2`.

This albedo is the most common used large-scale feedback in energy balance models, e.g.,
[Ghil1981](@cite), although it is typically taken as a piece-wise linear function.
There is little change with using a hyperbolic tangent instead, while the `tanh`
offers a differentiable flow.

The timescale `τ` if not zero will make an [`ExpRelaxation`](@ref) process
relaxing to the hyperbolic tangent.
"""
Base.@kwdef struct IceAlbedoFeedback <: Process
    max = 0.45
    min = 0.1
    Tscale = 10
    Tfreeze = C_to_K
    τ = 0 # if 0 you have no derivative
    α_ice = α_ice
    T = T
end

ProcessBasedModelling.lhs_variable(a::IceAlbedoFeedback) = a.α_ice
ProcessBasedModelling.timescale(a::IceAlbedoFeedback) = a.τ

function ProcessBasedModelling.rhs(a::IceAlbedoFeedback)
    y = a.α_ice
    suffixes = ["max", "min", "Tscale", "Tfreeze"]
    values = (a.max, a.min, a.Tscale, a.Tfreeze)
    mtk_pars = map((val, s) -> new_derived_named_parameter(y, val, s), values, suffixes)
    Tref = mtk_pars[4] - mtk_pars[3]/2
    # don't use `TanhProcess`, use the function instead
    iproc = ExpRelaxation(y, tanh_expression(a.T, mtk_pars[1:3]..., Tref), a.τ)
    return ProcessBasedModelling.rhs(iproc)
end
# we don't need to extend `lhs`, `ProcessBasedModelling` takes care of that.