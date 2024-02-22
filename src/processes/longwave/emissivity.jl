export EmissivityFeedbackTanh, EmissivitySellers1969, SoedergrenClearSkyEmissivity

"""
    EmissivityFeedbackTanh(; ε, Τ, left = 0.5, right = 0.4, rate = 0.5, Tref = 288.0)

Create an equation that assigns emissivity `ε` to hyperbolic tangent of temperature `T`.
This is an ad-hoc feedback  that was used in [Bastiaansen2023](@cite),
similar to [`EmissivitySellers1969`](@ref).
In essence this is a [`TanhProcess`](@ref) with the given keywords as parameters.
"""
function EmissivityFeedbackTanh(; left = 0.5, right = 0.4, rate = 0.5, Tref = 288.0)
    return TanhProcess(ε, T, left, right, rate, Tref)
end

"""
    SoedergrenClearSkyEmissivity(; ε=ε, T=T, RH = 0.8, H_H20 = 2.0)

Create Eq. 10 of [Soedergren2018](@cite), which is the same as Eq. 21 of [Barker1999](@cite)
for the effective emissivity of clear sky atmosphere:
```math
\\varepsilon = 1 - \\exp(0.082 - (2.38*0.1*e_s*RH*H_H2O + 40.3*CO2*1e-6)^0.294)
```
with ``e_s`` the [`saturation_vapor_pressure`](@ref).
The equation assumes CO2 concentration is in ppm and vapor pressure in kPa
hence the conversion factors 0.1 and 1e-6.

!!! warn "Physically wrong equation"
    Be advised: this process is included for reference only. It should not be used
    because it is physically wrong. Emissivity increases with temperature,
    while it should decrease: higher temperature → stronger greenhouse effect
    → smaller effective emissivity required for higher temperature to have
    the same OLR as per the basic equation ``OLR = ε σ Τ^4``.
"""
function SoedergrenClearSkyEmissivity(; ε = ε, T = T, RH = 0.8, H_H2O = 2.0)
    e_s = saturation_vapor_pressure(T)
    ε ~ 1 - exp(0.082 - (2.38*0.1*e_s*RH*H_H2O + 40.3*CO2*1e-6)^0.294)
end

"""
    EmissivitySellers1969(; ε = ε, T = T, m = 0.5)

Create the equation `ε ~ 1 - m*tanh(19*T^6*1e-16)`, which was used originally
in [Sellers1969](@cite) to represent the effect of "water vapor, carbon dioxide,
dust and clouds on terrestrial radiation".
"""
function EmissivitySellers1969(; ε = ε, T = T, m = 0.5)
    return ε ~ 1 - m*tanh(19*T^6*1e-16)
end
