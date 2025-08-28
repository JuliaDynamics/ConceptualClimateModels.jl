"""
    CloudToppedMixedLayerModel

Submodule providing processes about cloud topped mixed layer models (MLMs).
This combines existing equations on MLM by [Stevens2006](@cite) and [Bretherton1997](@cite),
with surface energy balance and dynamic cloud equations.
It is developed as part of the research article [Datseris2025](@cite).
If you use this submodule, please cite the paper.

The organization is as follows:

- All important variables and parameters (participate in many processes)
  are defined in the module file and as module-level scoped
  variables (global variables). Click the "source" of this docstring
  to access the module file.
- All processes (physical equations) are defined in their respective files
  such as `free_troposhere.jl`, etc. Docstrings of important processes
  are expanded in the docs here. You will notice that all functions that
  return processes (equations) utilize these global variables and global parameters.
  Many of these functions will also define local variables and parameters.
  All noteworthy processes have docstrings that are expanded in the submodule online
  documentation. However, the majority of docstrings do not actually list
  the equations themselves. Simply click the "source code" button on the bottom
  right of each docstring to go to the source code. Because this package is written
  in Julia, and because it uses symbolic expressions throughout, reading the source
  code is truly as straight forward as reading Latex-rendered equations.
- The `default.jl` file defines default processes for many global variables.
  These are also expanded in the docs.

Throughout the submodule time is in units of days, specific humidity is in units of g/kg,
liquid water static energy is in units of K (i.e., normalized by cₚ),
height in meters, temperature in K, and all energy quantities are in W/m².

To learn how to use this submodule visit first the general tutorial of ConceptualClimateModels.jl
and then the dedicated example on a cloudy mixed layer model.
The module purposefully does not export any names, so the recommended way to use it
is by an alias: `import ConceptualClimateModels.CloudToppedMixedLayerModel as CTMLM`.
"""
module CloudToppedMixedLayerModel

using ConceptualClimateModels
import Roots, NaNMath

###########################################################################################
# Constants, universal parameters, variables, and observables
###########################################################################################
const ℓ_v = 2.53e6 # latent heat of vaporization (J/gr) # notice the gram units! coz we use `q` in gr/kg!
const p₀ = 101780.0 # reference pressure at the sea surface in Pa
const e0 = 610.78 # reference water vapor partial pressure at 273.15 K in Pa = J/m^3
const Rd = 287.0 # gas constant dry air (J/K/kg)
const Rv = 461.0 # gas constant water vapor (J/K/kg)
const σ_SB = 5.670374419e-8 # stefan boltzman constant
const cₚ = 1004.0 # heat capacity at constant pressure (J/K/kg)
const g = 9.8 # gravity
const C_to_K = 273.15 # absolute zero temperature, add to Kelvin to make C
const sec_in_day = 60*60*24.0

# Universal parameters: environmental conditions. These are used in more than one functions
# that provide equations/processes/parameterizations.
ctbbl_parameters = @parameters begin # the default values give results very close to Stevens 2006
    # external conditions that could come from e.g., ERA5
    (U = 6.0), [bounds = (3.0, 9.0), description = "surface wind speed (m/s)"]
    # The drag coefficient is crucial. Smaller values promote more stratocumulus
    # while having a small direct impact on SST (through evaporative cooling).
    # Both `U` and `d_c` are equivalent in the eyes of the model; they have same effect.
    (d_c = 0.0012), [bounds = (7.9e-4, 0.0021), description = "aerodynamic bulk drag/transfer coefficient"]
    # values of V = d_c*U range from = [
    #     0.015, # Lilly 1968
    #     0.0021 * 7, # chung and teixeira 2012, same as Lilly
    #     8e-4*7.0 = 0.0056, # Singer & Schneider 2023a
    #     7.9e-4*10.0 = 0.0079, # Singer & Schneider 2023b
    #     0.0011*8 = 0.0088, # Stevens 2006
    #     0.0012*5.7 = 0.00684 # deRoode 2016
    #     While Bretherton and Wyant 1997 use d = 0.001V*(1 + 0.007V) (eq. 6) which is nonlinear in U
    #     and for U in (4, 8) it would give (0.004112, 0.008448)
    # ]
    (D = 3e-6), [description = "large scale divergence, 1/s"]
    (RH₊ = 0.2), [description = "relative humidity above boundary layer"]
    (δ_Δ₊T = 5.0), [description = "prescribed variability in temperature inversion (difference), K"]
    (δ_FTR = 0.0), [description = "Environmental variability of T_FTR, K"]
    (CO2 = 400.0), [description = "CO2 concentration, ppm"]
    (ECS_CO2 = 3.0), [description = "tropospheric warming per doubling of CO2 (global warming), K"]
    # timescales
    (τ_C = 2.0), [description = "exponential relaxation timescale of cloud fraction"]
    (τ_SST = 50.0), [description = "SST timescale in days (multiply with `sec_in_day` to get heat capacity in J/K/m²)"]
end

ctbbl_variables = @variables begin
    # Variables that can be dynamic (d/dt) have a default value.
    # Clouds
    (C(t) = 1.0), [bounds = (0.0, 1.0), description = "cloud fraction"]
    i_Λ(t), [description = "decoupling index: 0 at coupled, 1 at decoupled"]
    LWP(t), [description = "liquid water path within the cloud layer, g/m^2"]
    z_cb(t), [description = "height of the cloud base, m"]
    z_ct(t), [description = "height of the cloud top = top of the boundary layer = at inversion, m"]
    ε_C(t), [description = "emissivity of cloud layer, normally per cloud fraction"]
    α_C(t), [description = "albedo of cloud layer per cloud fraction"]
    T_C(t), [description = "cloud layer emission temperature, K"]
    L_c(t), [description = "Cloud-top emitted longwave radiation, W/m²"]
    C_Λ(t), [description = "cloud fraction (sigmoid) as a function of decoupling"]
    C_κ(t), [description = "thresholded cloud fraction as a function of CLT"]
    C_∞(t), [description = "relaxation cloud fraction"]
    Λ(t), [description = "decoupling pamarameter: when ≥ Λc, Sc decouples into Cu"]
    CLT(t), [description = "cloud layer thickness = (z_ct - z_cb)/z_b, m"]
    RCT(t), [description = "relative cloud thickness (normalized) = CLT/z_b"]
    CRC(t), [description = "cloud layer total radiative cooling, W/m²"]
    CTRC(t), [description = "cloud top total radiative cooling, W/m²"]
    CTRClw(t), [description = "longwave cloud-top radiative cooling, W/m²"]
    CRClw(t), [description = "longwave component of CRC, W/m²"]
    CRCsw(t), [description = "shortwave component of CRC, W/m²"]
    w_v(t), [description = "ventilation velocity in the absence of stratocumulus"]
    w_m(t), [description = "mass influx velocity in the absence of stratocumulus"]
    # Surface
    (SST(t) = 290.0), [bounds = (270.0, 310.0), description = "SCT sea surface temperature"]
    SST_X(t), [description = "Energy exchange at SST due to diffusion, ocean heat uptake, or whatever else"]
    LHF(t), [description = "surface latent heat flux, W/m²"]
    SHF(t), [description = "surface sensible heat flux, W/m²"]
    ASW(t), [description = "shortwave radiation absorbed by surface, W/m²"]
    α(t), [description = "total albedo"]
    Ld(t), [description = "downwards longwave radiation reaching the surface, W/m²"]
    L₀(t), [description = "upwards longwave radiation emitted by surface, W/m²"]
    Lnet(t), [description = "net longwave radiation leaving (cooling) the surface, W/m²"]
    L_b(t), [description = "down/up longwave radiation emitted by boundary layer, W/m²"]
    L_FTR(t), [description = "longwave radiation emitted by free troposphere, W/m²"]
    ρ₀(t), [description = "moist air density at the surface"]
    q₀(t), [description = "sea surface total water specific humidity = saturation specific humidity"]
    s₀(t), [description = "sea surface liquid water static energy"]
    # Boundary layer / thermodynamics
    (q_b(t) = 11.0), [bounds = (1.0, 30.0), description = "total water specific humidity of boundary layer, g/kg"]
    (s_b(t) = 290.0), [bounds = (280.0, 310.0), description = "liquid water static energy of boundary layer, K (cₚ units)"]
    (z_b(t) = 1200.0), [bounds = (600.0, 2000.0), description = "height of the top of the boundary layer, m"]
    RH_b(t), [description = "relative humidity of boundary layer (w.r.t. surface saturation)"]
    ΔF_s(t), [description = "energy forcing of the BL, W/m²"]
    ΔF_q(t), [description = "moisture forcing of the BL, W/m²"]
    T_t(t), [description = "temperature at the cloud top = top of the boundary layer = before inversion, K"]
    T_lcl(t), [description = "temperature at the lifting condensation level, K"]
    z_lcl(t), [description = "height of the lifting condensation level starting from surface, m"]
    w_e(t), [description = "entrainment velocity"]
    T_b(t), [description = "emission temperature of boundary layer"]
    Δ₊s(t), [description = "s₊ - s_b, as in Stevens2006"]
    Δ₊sᵥ(t), [description = "virtual s₊ - s_b, as in Singer Schneider"]
    Δ₊q(t), [description = "q₊ - q_b, as in Stevens2006"]
    Δ₀q(t), [description = "q_b - q₀, as in Stevens2006"]
    Δ₀s(t), [description = "s_b - s₀, as in Stevens2006"]
    ε_b(t), [description = "emissivity of bulk layer due to water vapor"]
    ε_t(t), [description = "total emissivity: ε_C + ε_β"]
    q_x(t), [description = "export of total humidity, g/kg/day"]
    s_x(t), [description = "export of static energy, K/day"]
    V(t), [description = "surface exchange velocity = U*drag"]
    ζ(t), [description = "cumulus state nondimensional mixing height as in Stevens 2006"]
    λ_s(t), [description = "decoupling ratio for liquid water static energy s, as in Eq.(1) of de Roode 2016"]
    λ_q(t), [description = "decoupling ratio for total water specific humidity"]
    # Above boundary layer/cloud top beyond
    T_FTR(t), [description = "emission temperature of free troposphere"]
    Δ₊T(t), [description = "temperature inversion strength, T₊ - T_t, K"]
    T₊(t), [description = "above boundary layer temperature (after inversion)"]
    s₊(t), [description = "above boundary layer liquid water static energy"]
    q₊(t), [description = "above cloud total water specific humidity"]
    ε_FTR(t), [description = "emissivity of the free troposphere"]
    S(t), [description = "insolation, W/m²"] # not a parameter, because can turn seasonal cycle on.
end

# include all other processes files
for (root, dirs, files) in walkdir(joinpath(@__DIR__))
    for file in files
        file == "CloudToppedMixedLayerModel.jl" && continue
        include(joinpath(root, file))
    end
end

end
