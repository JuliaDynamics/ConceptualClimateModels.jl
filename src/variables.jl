# These are quantities that change with time. Hence, they are either state variables
# or observables of the state variables. They have given names and are utilized
# in the rest of the files to make the equations that compose the dynamical system.
# In the final generated energy balance model it is not necessary that all of these
# will exist in the equations.

@variables T(t) = 290.0       # temperature, in Kelvin
@variables S(t) = 1.0         # insolation in units relative to solar constant
@variables f(t) = 0.0         # external forcing, normalized to units of the solar constant
@variables α(t) = 0.3         # planetary albedo, unitless
@variables α_ice(t) = 0.05    # albedo of ice, unitless
@variables α_cloud(t) = 0.1   # albedo of clouds, unitless
@variables ΔT(t) = 17.0       # equator to pole temperature difference, in Kelvin
@variables ε(t) = 0.5         # planetary effective emissivity, unitless
@variables C(t) = 0.6744      # cloud fraction, unitless
@variables ℓ(t) = 0.8         # sine of latitude of ice-line
@variables CO2(t) = 400       # CO2 concentration, in ppm

# Observables that can never be dynamic variables and hence have no default value:
@variables OLR(t) # outgoing longwave radiation
@variables ASR(t) # absorved shortwave radiation

export T, S, f, α, α_ice, α_cloud, ΔT, ε, ℓ, C, CO2, OLR, ASR

# This function is only meaningful for dynamic variables
function physically_plausible_limits(var::String)::Tuple{Float64, Float64}
    if var[1] == 'T'
        return (200, 350)
    elseif var == "α_ice" || var == "α_clouds"
        return (0, 0.75)
    elseif var[1] == 'α' || var[1] == 'ε' || var[1] == 'C'
        return (0, 1)
    elseif var == "ΔT"
        return (5.0, 60.0)
    else
        error("Unpsecified plausible physical limits for $(var). "*
        "Please edit function `physically_plausible_limits` and add one.")
    end
end
