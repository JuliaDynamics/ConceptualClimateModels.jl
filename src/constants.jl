const p_reference = 100000.0  # reference pressure in Pa
const triple_point = 273.16   # triple point (K)
const e0 = 610.78             # reference water vapor partial pressure (Pa = J/m^3)
const L0 = 2.5e6              # latent heat of vaporization (J/kg)
const Rv = 461.0              # gas constant water vapor (J/K/kg)
const Rd = 287.0              # gas constant dry air (J/K/kg)
const g = 9.8                 # gravity (m/s^2)
const cₚ = 1004.0              # heat capacity at constant pressure (J/K/kg)
const sec_in_day = 60*60*24.0
const C_to_K = 273.15 # absolute zero temperature, add to Kelvin to make C

const σ_Stefan_Boltzman = 5.670374419e-8    # stefan boltzman constant
const σ_Stefan_Boltzman_4root = 0.015431060170831203    # to the fourth root (for numerical stability)
const solar_constant = 340.25 # W/m^2, already divided by 4
const normalized_stefan_boltzmann_4root = (σ_Stefan_Boltzman/solar_constant)^0.25
const sec_in_year = 365*24*60*60.0

export sec_in_year, sec_in_day, C_to_K