"""
    cf_dynamic()

Provide the equations ``\\tau_C dC/dt = C_\\mathcal{D} - C`` as well as the
equation that defines ``C_\\mathcal{D}`` as a function of the decoupling
variable ``\\mathcal{D}``. The function uses the curve fitted to data
in [Datseris2025](@cite).
"""
function cf_dynamic(fit = :sigmoid)
    # During the reserch of the project I did a bunch of different fits.
    # The paper shows only the sigmoidal fit.

    starts = Dict(:exp => 1.0, :power => 1.0, :sigmoid => 0.5)
    scales = Dict(:exp => 0.4, :power => 0.8, :sigmoid => 1.0)

    @parameters begin
        (𝒟t = starts[fit]), [description = "scale of 𝒟 over which we go to Stratocumulus saturation in the sigmoidal curve"]
        (𝒟s = scales[fit]), [description = "when 𝒟>𝒟t the boundary layer is practically decoupled and we transition to cumulus"]
        (Cmax = 1.0), [description = "maximum (stratocumulus) cloud fraction"]
        (Cmin = 0.2), [description = "minimum (cumulus) cloud fraction"]
    end

    if fit == :exp
        function fitexp(D, Cmax, Cmin, start, scale)
            # Because of using MTK, it is better to clamp than to short circuit
            # D < start && return Cmax
            return min((Cmax - Cmin)*exp(-scale*(D - start)) + Cmin, Cmax)
        end
        C_𝒟_proc = C_𝒟 ~ fitexp(𝒟, Cmax, Cmin, 𝒟t, 𝒟s)
    elseif fit == :power
        C_𝒟_proc = C_𝒟 ~ ifelse(D < 𝒟t, Cmax, 1/(D-𝒟t+1)^𝒟s)
    elseif fit == :sigmoid
        C_𝒟_proc = SigmoidProcess(C_𝒟, 𝒟; left = Cmax, right = Cmin, scale = 𝒟s, start = 𝒟t)
    else
        error("unknown specification")
    end

    return [
        C_𝒟_proc,
        ExpRelaxation(C, C_𝒟, τ_C),
        # the decoupling index is just a translation of decoupling fit
        i_𝒟 ~ (Cmax - C_𝒟)/(Cmax - Cmin),
    ]
end

"""
    cloud_emission_temperature(version = :mean)

Return a process for ``T_C``. Versions are `:top, :base, :mean`.
"""
function cloud_emission_temperature(version = :mean)
    if version == :mean
        return T_c ~ (T_t + T_lcl)/2
    elseif version == :top
        return T_c ~ T_t
    elseif version == :base
        return T_c ~ T_lcl
    end
end

"""
    cloud_emissivity(version = 1.0; fraction = true)

Provide an equation for the effective emissivity of the cloud layer.
Options for `version`:

- `:clt`: inspired by [RandalSuarez1984](@cite), emissivity scales with the depth
  of the cloud layer.
- `:lwp`: Expression given by [Stephens1978](@cite) where emissivity is an exponential of LWP.
- `<: Number`: emissivity is just the provided number or symbolic expression.

If `fraction = true` the emissivity is further multiplied by the cloud fraction.
"""
function cloud_emissivity(version = 1.0; fraction = true)
    if version == :fraction
        expr = 1.0
    elseif version == :clt
        # inspired by Randal & Suarez 1984, emissivity depends on "depth"
        # The depth 200 is chosen to match the 10mb pressure used in Randal & Suarez
        @parameters ε_c_depth = 100.0 [description = "depth above which ε_c becomes 1"]
        expr =  min(CLT/ε_c_depth, 1) # use smoothstep if you don't want clamping
    elseif version == :lwp
        LWP = liquid_water_path_constql()
        expr = 1 - exp(-0.158*LWP)
    elseif version isa Number
        expr = version
    else
        error("incorrect version")
    end
    if fraction
        expr *= C
    end
    return ε_c ~ expr
end

function liquid_water_path_constql(T_t = T_t, z_cb = z_lcl, z_ct = z_b, q_b = q_b) # cloud base and top heights
    if z_cb ≥ z_ct
        return 0.0
    elseif any(isnan, (z_cb, z_ct))
        return NaN
    end
    # we have to do the integral of ρ*q_l over z from z_lcl to z_b.
    # To do so we will use the assumption that q_l increases linearly
    # from being 0 at the cloud base, to its max value at the cloud top.
    # Estimating the max value is rather simple:
    q_l_top = q_liquid(T_t, q_b, z_ct)
    # we will also use the assumption that the density remains constant
    # in the height of the cloud which allows us to analytically resolve the integral
    # (otherwise it is a function of temperature and the integral cannot be resolved)
    ρ_ref = moist_air_density(z_cb, T_t)/1e3 # correct units + use different height and temperature to get the average
    return 0.5*ρ_ref*q_l_top*(z_ct - z_cb)^2
end
@register_symbolic liquid_water_path_constql(T_t, z_cb, z_ct, q_b)

function liquid_water_path_exact(T_t, RCT, z_b, s_b, q_b)
    z_lcl = z_b - RCT*z_b # base of cloud layer
    if z_lcl ≥ z_b
        return 0.0
    elseif any(isnan, (z_b, RCT))
        return NaN
    end
    # normally we would use the exact temperature
    # T(z) = temperature_exact(z, s_b, q_b)
    # but it is safe to assume that temperature decreases linearly
    # within the cloud layer and with constant slope. We can estimate this
    # by obtaining the temperature at cloud top and cloud base,
    # assuming 0 liquid water at cloud base:
    T_b = s_b - g*z_lcl/cₚ # this is more than T_t
    # we then define the linear interpolation
    T(z) = Tb + (T_t - Tb)*(z - z_lcl)/(z_b - z_lcl)
    # We do exactly the same thing for the liquid water, which also increases linearly
    # and is zero at cloud base
    q_l_top = q_liquid(T_t, q_b, z_b)
    q_l(z) = q_l_top*(z - z_lcl)/(z_b - z_lcl)
    # We now do a discretized integral
    dz = 10.0
    zs = z_lcl:dz:z_b
    LWP = 0.0
    for z in zs
        LWP += moist_air_density(z, T(z)) * q_l(z) * dz
    end
    # Note that this version has convergence problems in the ODE solve
    # the resulting curve LWP(t) vs t is not smooth
    return LWP
end
@register_symbolic liquid_water_path_exact(T_t, RCT, z_b, s_b, q_b)


###########################################################################################
# Cloud base height / lifting condensation level
###########################################################################################
"""
    cloud_base_height(version = :exact, z_cb = z_lcl)

Provide an equation for the cloud base height captured by variable `z_cb`.
- `:exact`: exact estimation by figuring out when `q_liquid` first becomes positive.
  Computationally costly as it requires interpolations.
- `:Bolton1980`: Well known approximate expression by Bolton, 1980.

Because so far all versions calculate the lifting condensation level, `z_cb`
defaults to `z_lcl`. (and the default process for `z_cb` is for it to be `z_lcl`).
"""
function cloud_base_height(version = :exact, z_cb = z_lcl)
    if version == :exact
        z = cloud_base_height_exact(s_b, q_b, z_b)
    elseif version == :Zhang2006 # same as 2005 and 2009 papers and 2006 dissertation
        z = cloud_base_height_zhang2006(s_b, q_b, z_b)
    elseif version == :Bolton1980
        z = cloud_base_height_bolton1980(s_b, q_b)
    else
        error("incorrect specification for type for the lcl")
    end
    return z_cb ~ z
end

function cloud_base_height_exact(s, q, z)
    # Uses the algorithm described by Stevens 2006 when defining equation (39),
    # which for whatever reason does not take into account that pressure decreases with height.
    # First, assume q_liquid = 0 for finding the LCL.
    # This means that T = s - g*z/cₚ, since s is in units of cₚ
    f(x) = q - q_saturation(s - g*x/cₚ, x)
    if f(z) ≤ 0 # there is no liquid water
        z_lcl = z
    elseif f(0.0) ≥ 0 # very weird if this happens...
        z_lcl = 0.0
    else
        # Here we need to use a root finding algorithm.
        # # The algorithm `A42` appears to be the fastest option, even faster than bisection
        rootprob = Roots.ZeroProblem(f, (0.0, z))
        z_lcl = Roots.solve(rootprob, Roots.A42())
        # But the above sometimes leads to "no bracketing interval"
        # so I use instead
        # rootprob = Roots.ZeroProblem(f, z/2)
        # z_lcl = Roots.solve(rootprob)
    end
    return z_lcl
end
@register_symbolic cloud_base_height_exact(s, q, z)

# This does not give identical results to the "exact" version,
# but it is sliiightly higher. But really not by much!
"""
    cloud_base_height_bolton1980(s, q) → LCL (meters)

Compute the Lifting Condensation Level (LCL) for a given
liquid water static energy and specific humidity, from which
surface air temperature and relative humidity are extracted.

This is height (relative to parcel height) at which the parcel would become saturated
during adiabatic ascent. It is based on an approximate formula from Bolton
(1980 MWR) as given by Romps (2017 JAS).
For an exact formula see Romps (2017 JAS), doi:10.1175/JAS-D-17-0102.1
"""
function cloud_base_height_bolton1980(s, q)
    # given bulk conditions find surface temperature
    T = temperature_0_z(s)
    # find relative humidity at the surface
    RH = clamp(q / q_saturation(T), 0.0, 1.0)
    Tadj = T - 55.0  # in Kelvin
    return cₚ/g*(Tadj - (1/Tadj - log(RH)/2840.0)^(-1))
end
@register_symbolic cloud_base_height_bolton1980(s, q)


# This gives identical results to the Bolton version
# while being dramatically more expensive computationally.
# We do not use it anywhere therefore.
# (That's why `LambertW` is not installed, if you try to run it it will error)
function cloud_base_height_romps2017(s, q)
    T = temperature_0_z(s) # equal to s.
    p = p₀
    rh = clamp(q / q_saturation(T), 0.0, 1.0)
    # The rest comes from https://romps.berkeley.edu/papers/pubdata/2016/lcl/lcl.py
    # but keeping only the water branches
    Ttrip = 273.16     # K
    ptrip = 611.65     # Pa
    E0v   = 2.3740e6   # J/kg
    ggr   = 9.81       # m/s^2
    rgasa = 287.04     # J/kg/K
    rgasv = 461        # J/kg/K
    cva   = 719        # J/kg/K
    cvv   = 1418       # J/kg/K
    cvl   = 4119       # J/kg/K
    cpa   = cva + rgasa
    cpv   = cvv + rgasv
    function pvstarl(T)
        ptrip * (T/Ttrip)^((cpv-cvl)/rgasv) *
        exp( (E0v - (cvv-cvl)*Ttrip) / rgasv * (1/Ttrip - 1/T) )
    end

    pv = rh * pvstarl(T)
    qv = rgasa*pv / (rgasv*p + (rgasa-rgasv)*pv)
    rgasm = (1-qv)*rgasa + qv*rgasv
    cpm = (1-qv)*cpa + qv*cpv
    if rh == 0
       return cpm*T/ggr
    end
    aL = -(cpv-cvl)/rgasv + cpm/rgasm
    bL = -(E0v-(cvv-cvl)*Ttrip)/(rgasv*T)
    cL = pv/pvstarl(T)*exp(-(E0v-(cvv-cvl)*Ttrip)/(rgasv*T))
    lcl = cpm*T/ggr*( 1 - bL/(aL*real(LambertW.lambertw(bL/aL*cL^(1/aL), -1))) )
    return lcl
end
@register_symbolic cloud_base_height_romps2017(s, q)
