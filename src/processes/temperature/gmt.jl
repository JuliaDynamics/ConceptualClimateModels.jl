export BasicRadiationBalance
export heat_capacity_in_seconds

abstract type Temperature <: Process end
ProcessBasedModelling.lhs_variable(p::Temperature) = p.T
ProcessBasedModelling.timescale(proc::Temperature) = heat_capacity_in_seconds(proc.c_T)
heat_capacity_in_seconds(c_T) = c_T/solar_constant

"""
    BasicRadiationBalance(; T, f, c_T = 5e8)

Create the equation
```math
c_T \\frac{dT}{dt} = S*(1 - α) - OLR + f
```
representing the most basic radiative energy balance at the top of the atmosphere
setting a global mean temperature,
see e.g., any introductory article [North1981](@cite), [Ghil1981](@cite).
`c_T` is the heat capacity of the system in J/K/m². However, for convenience,
the parameter added to the final equation is `τ_T` which is the timescale in seconds,
i.e., `c_T/solar_constant`. `S` is the received insolation,
, by default equal to `solar_constant`, but could e.g., be any astronomical forcing
such as [`AstronomicalForcingDeSaedeleer`](@ref).
`α` is the albedo and `f` any additional forcing such as [`CO2Forcing`](@ref).

The observable `ASR ~ S*(1 - α)` for the absorbed solar radiation
is defined, while `OLR` needs to be itself parameterized
(by default [`EmissivityStefanBoltzmanOLR`](@ref) is used).
"""
Base.@kwdef struct BasicRadiationBalance <: Temperature
    c_T = 5e8   # units of J/K/m2
    f = f # units of solar constant!
    T = T
end

function ProcessBasedModelling.rhs(proc::BasicRadiationBalance)
    return (ASR - OLR + proc.f)/solar_constant # divide with solar due to timescale in seconds
end
