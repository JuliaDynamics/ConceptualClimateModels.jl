export ΔTLinearRelaxation

# Equation (2) from
# The latitudinal temperature gradient and its climate dependence as inferred from foraminiferal δ18O over the past 95 million years
# https://doi.org/10.1073/pnas.2111332119
# here defined as the tempereature difference between (0, 30) and
# (60, 90) latitudes. My own research showed that the point you average
# does not matter for linearity, only for coefficients.
reference_equator_to_pole(T, A = 36.53, B = 0.658) = -B*(T - C_to_K) + A # convert to celcius


"""
    ΔTLinearRelaxation(; ΔT, T, τ = 5e6, A = 36.53, B = 0.658)

Create the equation
```math
\\tau_{\\Delta T}\\frac{\\Delta T}{dt} = \\Delta T_{ref}(T) - \\Delta T
```
which exponentially relaxes the equator-to-pole temperature difference `ΔT` to
its reference value
``\\Delta T_{ref}(T) = A - B*(T - 275.15)``, i.e., it decreases linearly with global mean
temperature ``T`` (in Kelvin).
The default values for `A, B` are obtained from Equation (2) of [Gaskell2022](@cite).
We also fitted paleoclimate data of [Osman2021](@cite) and found very similar results,
`A = 35.8, B = -1.11` for north hemisphere and `A = 27.4, b = -0.513` for south.

Here `ΔT` is defined as the temperature difference between average temperatures at (0, 30)
and (60, 90) latitudes. The timescale is taken as 2 months, although if
`τ = 0` is given, the equation ``\\Delta T ~ \\Delta T_{ref}(T)`` is created instead.
"""
function ΔTLinearRelaxation(; ΔT, τ = 5e6, A = 36.53, B = 0.658)
    ExpRelaxation(ΔT,  -B*(T - C_to_K) + A, τ)
end
