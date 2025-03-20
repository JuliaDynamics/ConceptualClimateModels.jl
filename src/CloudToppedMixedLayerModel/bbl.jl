"""
    bbl_dynamic()

Provide equations 1-3 in [Datseris2025](@cite) (or, 31-33 in [Stevens2006](@cite))
defining the bulk boundary layer dynamics.
An additional auxilary velocity ``-w_m`` is added in the equation for ``z_b``
and two auxilary export terms ``q_x, s_x`` are added to the equations for ``q_b, s_b``.
All these auxilarity terms are 0 by default (otherwise, assign a process to them).
"""
function bbl_dynamic()
    return [
        # We add two more velocities that can depend on clouds: ventillation due to cloud
        # clearing, and mass influx, again due to cloud clearing.
        TimeDerivative(z_b, w_e - D*z_b - w_m, LiteralParameter(1/sec_in_day)),
        # Stevens is a bit sloppy with units here because a multiplication with a density
        # is required. Specifically, ds/dt does not have the same units as ΔF/z.
        # one needs to further divide ΔF/z with a density ρ in units of kg/m^3.
        # I believe Stevens 2006 assumes a density of 1 and doesn't put it anywhere.
        # Here we do modify the equations with a density of ρ₀ dividing ΔF.
        # (Also note that s is in cₚ units, q in g/kg)
        TimeDerivative(s_b, (w_e*Δ₊s + SHF/ρ₀/cₚ - ΔF/ρ₀/cₚ)/z_b - s_x/sec_in_day, LiteralParameter(1/sec_in_day)),
        TimeDerivative(q_b, (w_e*Δ₊q + LHF/(ρ₀*(ℓ_v/1e3)))/z_b - q_x/sec_in_day, LiteralParameter(1/sec_in_day)),
    ]
end

"Provide equations that prescribe boundary conditions, so that when combined with
`bbl_dynamic` we get results indentical to Sec. 4.1. of [Stevens2006](@cite)."
function bbl_boundary_stevens2006()
    @parameters begin
        ΔF_fixed = 40.0 # this is W/m^2.
        q₊_fixed = 1.56 # e-3
        s₊_fixed = 301200.0/cₚ
    end
    # Stevens 2006 provides values for q₊, s₊ in Fig. 1.
    # I prescribe q₀ by providing a fixed SST as in the paragraph above equation (9)
    # The rest of the values are in Sec. 4.2.
    # where every component of the equations becomes a fixed constant:
    return [
        w_e ~ e_e*ΔF/(Δ₊s*cₚ), # ΔF needs to be divided with ρ0 but stevens 2006 assumes ρ0 = 1.
        s₊ ~ s₊_fixed,
        s₀ ~ s₊ - 12.5e3/cₚ,
        q₊ ~ q₊_fixed,
        ΔF ~ ΔF_fixed,
        ρ₀ ~ 1, # this is assumed by Stevens 2006 and is required to give exactly same results
        q₀ ~ q_saturation(288.96 + 1.25), # see discussion above Eq. (9) in Stevens 2006
    ]
end

"""
    bbl_steadystate_stevens2006(fixed; z_b, q_b, s_b, CLT)

Return the equations 35-38 in [Stevens2006](@cite)
describing the analytically solved steady state of the MLM.
These equations could be coupled to other parts of module but
we have a problem of circular dependency for the steady state of ``z_b``.
If we attempt to couple them with the dynamic equations for ``C``,
then the following:
```
z_b ~ h⃰ * (e_e*σ_38)/(1 + σ_38 - e_e),
σ_38 ~ V*Δs*cₚ/(ΔF/ρ₀),
```
yields a circular dependency: `z_b` depends on `σ_38` which depends on `ΔF`
which depends on `T_t` which depends on `z_b`.
To resolve this a `fixed` option is given, which can be any of:
`ΔF, z_b, w_s, T_t`. This quantity is set fixed and becomes a parameter
so that the equation for `z_b` is closed.
"""
function bbl_stevens2006_steadystate(fixed = :ΔF; z_b=z_b, s_b=s_b, q_b=q_b, CLT=CLT) # changing the variables allows defining observables!
    @variables σ_38(t), Δs(t), Δq(t), h⃰(t), η⃰(t), η_b(t)
    eqss = [

        # Let's first write the remaining equations:
        s_b ~ s₀ - Δs*(1 - e_e)/σ_38,
        q_b ~ q₀ + (q₊ - q₀)*(e_e)/(1 + σ_38),
        Δs ~ s₊ - s₀,
        Δq ~ q₊ - q₀,
        σ_38 ~ V*Δs*cₚ/(ΔF/ρ₀),  # we added the extra division with ρ₀ for correct units
        h⃰ ~ (ΔF/ρ₀)/(D*Δs*cₚ),   # we added the extra division with ρ₀ for correct units
        η⃰ ~ Rv*s₀^2/(ℓ_v*cₚ*g), # should be around 1500 m.
        η_b ~ -η⃰*log(1 + e_e/(1 + σ_38)*Δq/q₀) - ΔF/(V*g)*(1 - e_e),
        CLT ~ 1 - η_b/z_b,
        # We also need equations/values for q₊ and s₊
        # but this is done by combining these equations with `ftr_bblm_coupler`!
    ]

    # To get rid of the circular dependency in z_b we can fix either T_t
    # and deduce z_b to a value corresponding to an average lapse rate,
    # or fix z_b to some indicative value like 1000.0 and obtain
    # `T_t` from its formula. Or, we can even fix a subsidence velocity!
    if fixed == :z_b
        push!(eqss, ParameterProcess(z_b, 1400.0))
    elseif fixed == :T_t
        @parameters Γ_b = 0.008
        push!(eqss, ParameterProcess(T_t, 280.0), z_b ~ Γ_b*(SST - T_t))
    elseif fixed == :w_s
        @parameters w_s = 0.005 # from 0.001 to 0.01
        push!(eqss, z_b ~ w_e - w_s)
    elseif fixed == :ΔF
        push!(eqss,
            z_b ~ h⃰ * (e_e*σ_38)/(1 + σ_38 - e_e),
            ParameterProcess(ΔF, 40.0)
        )
    end
    return eqss
end

"""
    temperature_exact(z, s, q)

Use root-finding to find the temperature at height `z`
given the liquid water static energy and total specific humidity,
as described by [Stevens2006](@cite).
This is the default equation used for `T_t ~ temperature_exact(z_b, s_b, q_b)`.
"""
function temperature_exact(z, s, q)
    # actual moist static energy:
    s_actual(T) = T + g*z/cₚ - (ℓ_v/1e3)*q_liquid(T, q, z)/cₚ # units of g/kg for `q` hence /1e3
    # Temperature is by construction when static energy is actual static energy:
    f(T) = s - s_actual(T)
    Tguess = s - g*z/cₚ
    # Notice that here we use the "problem" interface
    # of CommonSolve.jl, which defines returns NaN on
    # non-convergence, instead of erroring.
    # This allows "bad" initial conditions to be normally
    # evolved and the ODE solving to hault.
    rootprob = Roots.ZeroProblem(f, Tguess)
    T = Roots.solve(rootprob, Roots.Order1(); atol = 0.01)
    return T
end
@register_symbolic temperature_exact(z, s, q)


##########################################################################################
# Entrainment
##########################################################################################
"""
    entrainment_velocity(version = :Stevens2006; use_augmentation = true)

Return an equation for the entrainment velocity ``w_e``.
Versions are `:Stevens2006, :Gesso2014, :LL96`.
Keyword `use_augmentation` adds the decoupling-based augmentation described in [Datseris2025](@cite).
"""
function entrainment_velocity(version = :Stevens2006;
        e_e = version == :Stevens2006 ? 0.9 : 0.5,
        w_e = w_e, ΔF = ΔF, Δs = Δ₊s, use_shear = false, use_augmentation = true
    )
    @parameters (e_e = e_e), [description = "Entrainment efficiency. Scales entrainment velocity"]

    if version == :Stevens2006
        x = e_e*ΔF/ρ₀/(Δs*cₚ) # corrected Stevens2006 with ΔF divided with ρ0
    elseif version ∈ (:LL96, :Gesso2014) # From Dal Gesso 2014
        # we obtain the beta coefficiencs from Stevens 2002 and set them to constants
        ϵ1 = 0.608
        Ad = βsd = 1 + ϵ1*q_b
        Bd = βqd = ϵ1
        As = βsm = 0.5
        Bs = βqm = 3.5
        # we then obtain the entrainment from Dal Gesso 2014 using the main
        # formula of the appendix for Lewellen and Lewellen (1998).
        # The factors of the equation are the coeffients after eq. (14) which
        # are based on the beta parameters (stated as A, B after equation 9).
        SLT = 1 - CLT # = cloud base height / inversion height
        # liquid water potential temperature inversion jump can be
        # estimated analytically as a function of s, q: θl = (s - g*z)/(cₚ*Π(z))
        # since the inversion jump height is approximated as 0,
        # Δ₊θl = (Δ₊s - 0)/(cₚ*Π(z_b)); and since our units of s are already
        # divided by cₚ, we have the very simple relationship of
        p_z = pressure(z_b, (T_t + T₊)/2) # so that I don't deal with change in pressure across inversion
        Π = (p_z/p₀)^0.286 # approximate Rd/cp; checked for correctness
        Δ₊θl = Δ₊s/Π
        Δ₊θvd = Ad*Δ₊θl + Bd*Δ₊q
        Δ₊θvs = As*Δ₊θl + Bs*Δ₊q
        if version == :LL96
            # This is the Lewellen & Lewellen 1998 version:
            x = 2*e_e*ΔF / (SLT^2*Δ₊θvd + (1 - SLT^2)*Δ₊θvs)
        elseif version == :Gesso2014
            # And this is equation 15 of Dal Gesso 2014:
            x = 5*e_e*ΔF / (2Δ₊θl + 2.5*e_e*(SLT^2*Δ₊θvd + (1 - SLT^2)*Δ₊θvs))
        end
        # in both cases the units of `x` appear to be already in mm/s so
        x = x/1e3
        # Unfortunately, during integration this equation for x is very unstable.
        # it can become negative due to the fact that Δ₊q is negative and can be
        # very negative as integration occurs. I made the arbitrary choice here
        # to bound x within reasonable limits to make integration more stable
        x = clamp(x, 0.0, 0.01) # Dal Gesso 2014 finds x up to 0.004 only.
    else
        error("incorrect entrainment")
    end
    if use_shear
        shear = 0.0005*exp(0.2)*exp(-z_b/500)
    else
        shear = 0
    end
    w = x + shear
    if use_augmentation
        @parameters β₋ = 0.3
        w = w*(1 + i_𝒟*CLT*β₋/2*V)
    end
    return w_e ~ w
end

##########################################################################################
# Thermodynamics
##########################################################################################
"""
    q_liquid(T, q, z = 0)

liquid specific humidity given temperature and total water specific humidity.
0.0 if below saturation.
"""
q_liquid(T, q, z = 0) = max(q - q_saturation(T, z), 0)

"When requesting temperature at height 0, we realize that
the expression for `s` is much simplified (no height, no liquid water)
so we have the temperature to be just `s` since we measure `s` in Kelvin."
temperature_0_z(s, args...) = temperature_0_z(s)
temperature_0_z(s) = s

"""
    q_saturation(T, z)

Saturation specific humidity given temperature and height.
Height is transformed to pressure via hydrostatic approximation.
"""
function q_saturation(T, z = 0)
    # Equation 39 of [Stevens2006](@cite).
    # Stevens never gives a value for q₀ so that we can use Eq. 39.
    # I used the values from: http://hvac-calculator.com/hum_pr2x.php
    # assumming 100% relative humidity, i.e., saturation.
    # Version from Bjorn Stevens 2006; There is something odd with this version
    # because it does not utilize height (pressure)...
    # Tref = C_to_K + 15
    # qref = 10.6
    # return qref*exp((ℓ_v/(Rv*Tref^2))*(T - Tref))
    # Version from Claire Singer codebase
    # uses Clasius-Clapeyron relation with assumed constant
    # latent heat of vaporization term L0=2.5e6
    psat = e0 * exp(-ℓ_v/Rv * (1/T - 1/273.16))
    return Rd/Rv * psat / (pressure(z, T) - psat)*1e3 # multiply with 1e3 due to my chosen units!

    # Version from Yunyan Zhang codebase (approximately same results)
    # p = pressure(z, T)
    # c0= 0.611213476e+03; c1= 0.444007856e+02; c2= 0.143064234e+01
    # c3= 0.264461437e-01; c4= 0.305930558e-03; c5= 0.196237241e-05
    # c6= 0.892344772e-08; c7=-0.373208410e-10; c8= 0.209339997e-13
    # eps = Rd/Rv
    # x = T-273.15 # to Celcius
    # # This is Horner's rule but alas we can't use the macro due to symbolic numbers :(
    # esl=c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
    # return eps*esl/(p-esl)*1e3 # multiply with 1e3 due to my chosen units!
end

"""
    pressure(z, T)

Use hydrostatic balance and ideal gas law to get pressure at height z
given temperature at height z. The equation is often called the "Hypsometric equation"
with the factor (Rd*T/g) called the scale height.
Note that normally using Rd requires usage of Tv (virtual temperature), defined as
Tv = T*(1 + 0.608*q) with q the specific humidity of water vapor.
In the codebase we practically always approximate Tv by T.
"""
pressure(z, T) = p₀ * exp((-g * z) / (Rd * T))

"moist air density"
function rho(z, T)
    return pressure(z, T) / (Rd * T)
end

"""
    temperature_constant_lapse(z, SST, Γ = 0.008)

Return the temperature at height `z` given `SST` and lapse rate `Γ`.
Uses the fact that lapse rate is approximately constant
and is between dry and moist, at around 8 - 8.5 K / km
(see e.g., Fig 5 of [Nowak2021](@cite)).
"""
function temperature_constant_lapse(z, SST, Γ = 0.008)
    # Clare Singer said that during the model evolution, before arrival at equilibrium
    # this is probably not an accurate thing to use
    return SST - Γ*z
end


"""
    moist_adiabatic_lapse_rate(T, z, RH)

Return the moist adiabatic lapse rate expected at given temperature,
height (to obtain pressure), and relative humidity.
"""
function moist_adiabatic_lapse_rate(T, z, RH)
    # get mixing ratio from https://github.com/Unidata/MetPy/blob/v1.6.2/src/metpy/calc/thermo.py#L2048-L2108
    epsilon = 0.622
    p_z = pressure(z, T)
    p_saturation = e0 * exp(-(ℓ_v/Rv) * (1/T - 1/C_to_K))
    w_s = saturation_mixing_ration = epsilon*p_saturation / (p_z - p_saturation)
    r = mixing_ratio = epsilon * w_s * RH / (epsilon + w_s*(1 - RH))
    # then get Lapse rate from
    # https://en.wikipedia.org/wiki/Lapse_rate#Moist_adiabatic_lapse_rate
    H_v = 2501000 # heat of vaporization of water =  J/kg
    R_sd =  287 # specific gas constant of dry air = J/kg·K
    R_sw = 461.5 # specific gas constant of water vapour = J/kg·K
    c_pd = 1003.5 #	the specific heat of dry air at constant pressure, =  J/kg·K
    Γ = g*(1 + (H_v*r)/(R_sd*T))/(c_pd + (H_v^2*r)/(R_sw*T^2))
    return Γ
end

