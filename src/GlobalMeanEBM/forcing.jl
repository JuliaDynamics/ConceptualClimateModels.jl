export CO2Forcing

"""
    CO2Forcing(; f, CO2, CO2f = 3.7)

Create the equation ``f ~ CO2f \\log_2(CO2/400)`` which describes the
forcing added to the TOA energy balance due to CO2 concentrations,
assumming the `OLR` expression is calibrated for 0 added forcing at 400 ppm
which is the default for OLR expressions provided by ConceptualClimateModels.jl.

The default value of ``f`` comes from Eq. (3.2) of [Bastiaansen2023](@cite) which cites IPCC-5,
while [Etminan2016](@cite) report practically the same value assuming a constant
``f`` (note here the log is base 2). In reality ``f`` depends on ``CO2`` and other
greenhouse gases concentrations due to spectral overlaps, see [Etminan2016](@cite) Sec. 4.
"""
Base.@kwdef struct CO2Forcing <: Process
    f = f
    CO2 = CO2
    CO2f = 3.7
end

ProcessBasedModelling.lhs_variable(::CO2Forcing) = f
ProcessBasedModelling.rhs(c::CO2Forcing) = c.CO2f*log2(CO2/400)