# These are quantities that change with time. Hence, they are either state variables
# or observables of the state variables. They have given names and are utilized
# in the rest of the files to make the equations that compose the dynamical system.
# In the final generated energy balance model it is not necessary that all of these
# will exist in the equations.

const GLOBAL_EBM_VARIABLES = @variables begin
    T(t) = 290.0       [description = "temperature, in Kelvin"]
    S(t) = 1.0         [description = "insolation in units relative to solar constant"]
    f(t) = 0.0         [description = "external forcing, normalized to units of the solar constant"]
    α(t) = 0.3         [description = "planetary albedo, unitless"]
    α_ice(t) = 0.05    [description = "albedo of ice, unitless"]
    α_cloud(t) = 0.1   [description = "albedo of clouds, unitless"]
    ΔT(t) = 17.0       [description = "equator to pole temperature difference, in Kelvin"]
    ε(t) = 0.5         [description = "planetary effective emissivity, unitless"]
    C(t) = 0.6744      [description = "cloud fraction, unitless"]
    ℓ(t) = 0.8         [description = "sine of latitude of ice-line"]
    CO2(t) = 400       [description = "CO2 concentration, in ppm"]
    # Observables that can never be dynamic variables and hence have no default value:
    OLR(t) [description = "outgoing longwave radiation"]
    ASR(t) [description = "absorved shortwave radiation"]
end

# TODO: Should we do this export...?
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
