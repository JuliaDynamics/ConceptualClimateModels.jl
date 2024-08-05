export EmissivityFeedbackTanh, EmissivitySellers1969, SoedergrenClearSkyEmissivity
export EmissivityLogConcentration

"""
    EmissivityFeedbackTanh(; ε, Τ, left = 0.5, right = 0.4, rate = 0.5, Tref = 288.0)

Create an equation that assigns emissivity `ε` to hyperbolic tangent of temperature `T`.
This is an ad-hoc feedback  that was used in [Bastiaansen2023](@cite),
similar to [`EmissivitySellers1969`](@ref).
In essence this is a [`TanhProcess`](@ref) with the given keywords as parameters.
"""
function EmissivityFeedbackTanh(; ε=ε, T=T, left = 0.5, right = 0.4, rate = 0.5, Tref = 288.0)
    return SigmoidProcess(ε, T, left, right, rate, Tref)
end

"""
    SoedergrenClearSkyEmissivity(; ε, T, CO2, RH = 0.8, H_H20 = 2.0)

Create Eq. 10 of [Soedergren2018](@cite), which is the same as Eq. 21 of [Barker1999](@cite)
for the effective emissivity of clear sky atmosphere:
```math
\\varepsilon = 1 - \\exp(0.082 - (2.38*0.1*e_s*RH*H_{H2O} + 40.3*CO2*1e-6)^{0.294})
```
with ``e_s`` the [`saturation_vapor_pressure`](@ref).
The equation assumes CO2 concentration is in ppm and vapor pressure in kPa
hence the conversion factors 0.1 and 1e-6.

!!! warn "Atmospheric, not effective emissivity!"
    Be advised: "effective emissivity" is a number multiplying surface outgoing
    radiation in models with only one layer, the surface: ε*σ*Τ^4. That is what
    other emissivity processes like [`EmissivitySellers1969`](@ref) represent.
    Here this atmospheric emissivity is actual, not effective. It is supposed to be included
    in two layer models with one layer being atmosphere and one being surface.
    That is why this ε here _increases_ with CO2/H2O concentrations.
    In contrast, the _effective_ emissivity of a surface-only model would
    decrease with CO2/H2O concentrations.
"""
function SoedergrenClearSkyEmissivity(; ε = ε, T = T, CO2 = CO2, RH = 0.8, H_H2O = 2.0)
    e_s = saturation_vapor_pressure(T)
    ε ~ 1 - exp(0.082 - (2.38*0.1*e_s*RH*H_H2O + 40.3*CO2*1e-6)^0.294)
end

"""
    EmissivitySellers1969(; ε, T, m = 0.5)

Create the equation `ε ~ 1 - m*tanh(19*T^6*1e-16)`, which was used originally
in [Sellers1969](@cite) to represent the effect of "water vapor, carbon dioxide,
dust and clouds on terrestrial radiation".
"""
function EmissivitySellers1969(; ε = ε, T = T, m = 0.5)
    return ε ~ 1 - m*tanh(19*T^6*1e-16)
end


"""
    EmissivityLogConcentration(; ε, T, CO2, q,
        ε_ref = 0.6, ε_CO2 = 0.1, ε_H20 = 0.01
    )

Create the equation
```math
ε ~ ε_{ref} - ε_{CO2}*log2(CO2/400) - ε_{H2O}*log2(q/1e-3)
```
This equation is inspired by Eq. A.9 of [VonderHeydt2016](@cite), which had only the
first two terms. Here we added the last term to represent feedback due to water
vapor. `q` is absolute humidity (or, water vapor concentration),
whose default process is [`AbsoluteHumidityIGLCRH`](@ref).

The logarithm of either water vapor or CO2 concentration is used because literature suggests
that the radiative forcing of greenhouse gases scales with the logarithm of their
concentration [Huang2014](@cite).
As `ε` is a first order component of the OLR, it needs to scale linearly with the logarithm.
"""
function EmissivityLogConcentration(;
        ε = ε, T = T, CO2 = CO2, q = q, # global variables used by default
        ε_ref = 0.6, ε_CO2 = 0.1, ε_H2O = 0.01, # with α=0.3 this gives temperature of today
    )
    # Water vapor contribution is similar to CO2, it comes from its concentration.
    # Absolute humidity is water concentration (up to multiplicative constant).
    @convert_to_parameters ε_ref ε_CO2 ε_H2O
    return ε ~ ε_ref - ε_CO2*log2(CO2/400.0) - ε_H2O*NaNMath.log2(q/1e-3)
end
