export saturation_vapor_pressure

"""
    saturation_vapor_pressure(T)

Given surface temperature (in K) return saturation pressure for water vapor (in kPa)
using the Tetens-Murray formula from [TetensEquation2023](@cite), which is `A*exp(B*T/(C+T))`.
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