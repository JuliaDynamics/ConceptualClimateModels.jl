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

"""
    ΔTStommelModel(; ΔT=ΔT, ΔS=ΔS, η1 = 2, η2 = 1, η3 = 0.3)

Create the equations
```math
\\dot{\\Delta T} = \\eta_1 - \\Delta T - |\\Delta T - \\Delta S| \\Delta T
\\dot{\\Delta S} = \\eta_2 - \\eta_3\\Delta S - |\\Delta T - \\Delta S| \\Delta S
```
which are the two equations of the Stommel box model for Atlantic thermohaline circulation
[Stommel1961](@cite), here presented in nondimensionalized form [Lohmann2021](@cite), so that
temperature and sality are normalized by their coefficients ``a_T, a_S`` relating them
to the density of water
```math
\\rho = \\rho_0 [1 - a_T(T - T_0) + a_S(S-S_0)]
```
for some reference values.
"""
function ΔTStommelModel(; ΔT=ΔT, ΔS=ΔS, η1 = 2, η2 = 1, η3 = 0.3)
    return [
        TimeDerivative(ΔT, η1 - ΔT - abs(ΔT - ΔS)*ΔT),
        TimeDerivative(ΔS, η2 - η3*ΔS - abs(ΔT - ΔS)*ΔS),
    ]
end