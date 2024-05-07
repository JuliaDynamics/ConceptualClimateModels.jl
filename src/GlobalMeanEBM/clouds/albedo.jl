export CloudAlbedoExponential, CloudAlbedoLinear

"""
    CloudAlbedoExponential(
        α_cloud, C, a = 2.499219232848238, b = 17.596369331717433
    )

Create the equation `α_cloud ~ sinh(a*C)/b` relating cloud albedo
to cloud fraction `C`. This equation is exponential and not linear, as in observations.
[Engstrom2015](@cite) (and also [Bender2017](@cite)) discuss this exponential relation
in detail, and provide as explanation that cloud effective albedo increases with latitude
(due to solar zenith changes) while cloud fraction also increases with latitude.

Note that here however we modify the equation `α_cloud ~ exp(a*C - b)` of [Engstrom2015](@cite)
to utilize the hyperbolic sine, so that `α_cloud = 0` when `C = 0` as is physically
necessary. Then, `a, b` are extracted by fitting CERES data, using as `α_cloud` the
energetically consistent cloud albedo as defined by [Datseris2021](@cite),
further yearly averaged and within latitudes (-60, 60) as in [Bender2017](@cite).
This albedo can be directly added to the clear sky albedo to produce the planetary albedo.
"""
function CloudAlbedoExponential(;
        α_cloud = α_cloud, C = C, a = 2.499219232848238, b = 17.596369331717433
    )
    return α_cloud ~ sinh(a*C)/b
end

"""
    CloudAlbedoLinear(; α_cloud, C, a = 0.2252861764703435)

Same as in [`CloudAlbedoExponential`](@ref), but now the linear form `α_cloud ~ a*C`
is returned, with `a` fitted from CERES data in the same way.
"""
function CloudAlbedoLinear(;
        α_cloud = α_cloud, C = C, a = 0.2252861764703435
    )
    return α_cloud ~ a*C
end
