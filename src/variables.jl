# module CCMV

# These are quantities that change with time. Hence, they are either state variables
# or observables of the state variables. They have given names and are utilized
# in the rest of the files to make the equations that compose the dynamical system.
# In the final generated energy balance model it is not necessary that all of these
# will exist in the equations.

const PREDEFINED_CCM_VARIABLES = @variables begin
    (T(t) = 290.0),      [bounds = (200.0, 350.0), description = "temperature, in Kelvin"]
    (S(t) = 1.0),        [bounds = (0.8, 1.2), description = "insolation, normalized to units of the solar constant"]
    (f(t) = 0.0),        [bounds = (-0.1, 0.1), description = "external forcing, in W/m²"]
    (α(t) = 0.3),        [bounds = (0.0, 1.0), description = "planetary albedo, unitless"]
    (α_ice(t) = 0.05),   [bounds = (0.0, 1.0), description = "albedo of ice and/or snow, unitless"]
    (α_cloud(t) = 0.1),  [bounds = (0.0, 1.0), description = "albedo of clouds, unitless"]
    (ΔT(t) = 17.0),      [bounds = (0.0, 60.0), description = "equator to pole temperature difference, in Kelvin"]
    (ΔS(t) = 5.0),       [bounds = (0.0, 10.0), description = "equator to pole salinity difference, in psu"]
    (ε(t) = 0.5),        [bounds = (0.0, 1.0), description = "planetary effective emissivity, unitless"]
    (C(t) = 0.6744),     [bounds = (0.0, 1.0), description = "cloud fraction, unitless"]
    (ℓ(t) = 0.8),        [bounds = (0.0, 1.0), description = "sine of latitude of ice-line"]
    (CO2(t) = 400.0),    [bounds = (200.0, 1800.0), description = "CO2 concentration, in ppm"]
    (q(t) = 1.0),        [bounds = (1e-5, 1e2), description = "absolute humidity, in gr/m³"]
    # Observables that can never be dynamic variables and hence have no default value:
    (OLR(t)), [description = "outgoing longwave radiation, in W/m²"]
    (ASR(t)), [description = "absorved shortwave radiation, in W/m²"]
end

export PREDEFINED_CCM_VARIABLES
export T, S, f, α, α_ice, α_cloud, ΔT, ΔS, ε, ℓ, C, CO2, OLR, ASR, q

function physically_plausible_limits(var::String)::Tuple{Float64, Float64}
    if var[1] == 'T'
        return (200, 350)
    elseif var[1] == 'α' || var[1] == 'ε' || var[1] == 'C'
        return (0, 1)
    else
        error("""
        Unpsecified plausible physical limits for $(var): it has no defined bounds or
        a default variable. You need to redefine the variable to have either of the two.
        """)
    end
end

# end