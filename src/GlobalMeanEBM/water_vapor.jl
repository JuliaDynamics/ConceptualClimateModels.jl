export saturation_vapor_pressure, AbsoluteHumidityIGLCRH, absolute_humidity_igl

"""
    AbsoluteHumidityIGLCRH(; q = q, T = T, RH = 0.5)

Create the equation
```math
q = m_v \\cdot RH \\cdot e_s / (R\\cdot T)
```
which approximates absolute humidity (water vapor concentration) in gr/m³
under two assumptions:
1. ideal gas law
2. constant relative humidity `RH` with warming, i.e., `RH` is a constant,
   so that `q = q(T)` _only_.

Here `m_v = 18.016` gr is the specific mass of water, `R` the ideal gas constant,
`e_s = `[`saturation_vapor_pressure`](@ref)`(T)` in kPa.

Note that you could create a process for a _variable_ `RH` and pass the variable
to the `RH` keyword to skip assumption (2).

Use `absolute_humidity_igl(T, RH)` to obtain `q` at a given `T, RH`.
"""
function AbsoluteHumidityIGLCRH(; q = q, T = T, RH = 0.5)
    @convert_to_parameters RH
    return q ~ absolute_humidity_igl(T, RH)
end

function absolute_humidity_igl(T, RH)
    f = mass_vapor * 1e3 / R_ideal_gas # sat vap pressure is in kPa, not Pa
    e_s = saturation_vapor_pressure(T)
    return f * e_s * RH / T
end

"""
    saturation_vapor_pressure(T)

Given surface temperature (in K) return saturation pressure for water vapor (in kPa)
using the Tetens-Murray formula from [TetensEquation2023](@cite), which is `A*exp(B*T/(C+T))`.
Different formula is used for when `T` is less than the freezing point.
"""
function saturation_vapor_pressure(T)
    T = T - C_to_K # the formula requires T in Celcius
    expr1 = let
        A = 0.61078
        B = 17.625
        C = 243.04
        A*exp(B*T/(C+T))
    end
    expr2 = let
        A = 0.61121
        B = 22.587
        C = 273.86
        A*exp(B*T/(C+T))
    end
    # This is done because we need to use `ifelse` for conditional
    # formulas within ModelingToolkit.jl.
    return ifelse(T>0, expr1, expr2)
end
