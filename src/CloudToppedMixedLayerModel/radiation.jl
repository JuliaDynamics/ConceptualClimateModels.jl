# All processes related to radiation: longwave emission, emissivities, albedo
# and cooling of the BL
# We make the assumption that the emission temperature of all layers is the same
# regardless of whether the radiation is emitted upwards or downwards.
# Especially inaccurate for the cloud layer as we use the cloud-top radiation for it!

#########################################################################################
# Clouds and BBL
#########################################################################################
"""
    cloud_shortwave_warming([version]; cloud_fraction = true)

Provide an equation for `CRCsw` (same as `CTRCsw`) which by default is `0.04*C*S`.
Otherwise `version` can be a `Number` specifying the RHS.
"""
function cloud_shortwave_warming(version = :insolation; cloud_fraction = true)
    CRCswrhs = if version isa Number
        version
    elseif version == :insolation
        0.04*S # 12 for S = 400, C = 0.75, which is the average in Zheng 2021
    elseif version == :constant
        12.0 # the evarge in Zheng 2021
    else
        error("incorrect specification for shortwave warming")
    end
    if cloud_fraction
        CRCswrhs *= C
    end
    return CRCsw ~ CRCswrhs
end

"""
    cloud_longwave_cooling(cloud_fraction = false)

Provide processes for `CTRClw, CRClw` based on the three-layer radiation balance.
If `cloud_fraction == false`, the RHS are further scaled by `C`.
By default this is not done because typically `ε_C` is scaled by `C`.
"""
function cloud_longwave_cooling(cloud_fraction = false)
    eqs = [
        CTRClw ~ L_c - ε_C*L_FTR, # cloud top cooling
        CRClw ~ CTRClw + L_c - ε_C*L_b - ε_C*(1 - ε_b)*L₀, # plus cloud bottom cooling
    ]
    if cloud_fraction
        eqs = map(eq -> eq.lhs ~ C*eq.rhs, eqs)
    end
    return eqs
end

"""
    mlm_radiative_cooling(version = :three_layer)

Provide an equation for ``\\Delta F_s``, the radiative cooling of the boundary layer
(assumming ``\\Delta F_q = 0``). Versions are: `:three_layer, :ctrc, :Gesso2014`
as in [Datseris2025](@cite).
"""
function mlm_radiative_cooling(version = :three_layer)
    if version == :ctrc
        return ΔF_s ~ CTRC
    elseif version == :crc
        return ΔF_s ~ CRC
    elseif version == :three_layer
        # I'll write stuff explicitly as I've fucked up before
        L_up_surf = L₀
        L_up_lcl = L_up_surf*(1 - ε_b) + L_b
        L_up_top = L_up_lcl*(1 - ε_C) + L_c

        L_down_top = L_FTR
        L_down_lcl = L_down_top*(1 - ε_C) + L_c
        L_down_surf = L_down_lcl*(1 - ε_b) + L_b

        return ΔF_s ~ (L_up_top - L_down_top) + (L_down_surf - L_up_surf)
    elseif version == :top_bottom
        bottom = L_b - ε_b*L₀
        ΔF_s ~ CTRC + bottom
    elseif version == :Singer2023
        # The equations in the paper are different
        # from the equations in the code... Plus it is impossible for the Singer
        # ΔTₑ to become -25 with this equation, as it can be at most -10.1
        # because q₊ is realistically never below 1.
        # So how do they plot -25 in the paper figure????????
        ΔF_s ~ C*0.9*σ_SB*(T_C^4 - (T_C - 10.1 + 3.1*log(CO2) + 5.3*log(q₊))^4)
    elseif version == :Gesso2014
        ΔF_s ~ max(82.0 - 7.9*q₊, 1.0) # I am clamping here for numerical stability. Perhaps I shouldn't?
    else
        error("Incorrect version for ΔF_s.")
    end
end

#########################################################################################
# Surface
#########################################################################################
"""
    downwards_longwave_radiation([version])

Provide equation for ``L_d, L_{net}``, the incoming longwave radiation or the net
longwave radiative cooling of the surface.
See the source code for possible versions.
"""
function downwards_longwave_radiation(version = :three_layer)
    if version == :three_layer
        return Ld ~ L_b + (1 - ε_b)*L_c + (1 - ε_b)*(1 - ε_C)*L_FTR
    elseif version == :two_layer
        return Ld ~ L_b + (1 - ε_b)*L_c
    elseif version == :single_layer
        return Ld ~ L_b
    elseif version == :fixed
        @parameters Lnet_0 = 56.0 [description = "Net LW cooling at surface"]
        return [Lnet ~ Lnet_0, Ld ~ L₀ - Lnet]
    elseif version == :clouds_only
        @parameters Lnet_b = 60 [description = "Fixed LW cooling due to BBL at surface"]
        return [Lnet ~ Lnet_b - (1 - ε_b)*L_c, Ld ~ L₀ - Lnet]
    else
        error("incorrect version for downwards longwave.")
    end
end

function bbl_emission_temperature(version = :parameterized)
    if version == :parameterized
        @parameters h_b = 0.5 [description = "Relative height of LW emission of BBL"]
        T_b ~ h_b*T_lcl + (1 - h_b)*s_b # s_b == temperature at surface
    elseif version == :surface
        T_b ~ s_b
    else
        error("unknown version")
    end
end

#########################################################################################
# Emissivity and albedo
#########################################################################################
"""
    matsunobu_emissivity(RH, T, H = 0)

Return emissivity at given relative humidity and temperature as defined by [Matsunobu2024](@cite).
"""
function matsunobu_emissivity(RH, T, H = 0) # default H means don't use height
    p = normalized_vapor_pressure(RH, T)
    h = 0.15*(exp(-H/8504) - 1) # height correction, by default zero
    ε = clamp(0.6 + 1.642*NaNMath.sqrt(p) + h + 0.02, 0, 1)
    return ε
end
function mendoza_emissivity(RH, T) # an alternative version I've considered but ignored
    p = normalized_vapor_pressure(RH, T)
    ε = clamp(1.108*p^(0.083), 0, 1)
    return ε
end

function logarithmic_emissivity(q) # somewhat ad hoc, I pulled this out of a hat. We know emissivity ∝ log(concentration)
    @parameters begin
        (ε_a = 0.7), [description = "reference emissivity of atmosphere at 1g/kg H2O specific humidity, and no clouds."]
        (ε_q = 0.05), [description = "increase in atmospheric emissivity due to water vapor (multiplies log(q))"]
    end
    return clamp(ε_a + ε_q*NaNMath.log(q), 0, 1)
end

function normalized_vapor_pressure(RH, T)
    # Use same formula as in the Matsunobu paper
    p_saturation = 610.94*exp(17.625(T − 273.15)/(T − 30.11))
    Pw = RH*p_saturation
    pw = Pw/p₀
    return pw
end

"Return an equation for `α` the total albedo perceived by the surface."
function albedo(version = :multiplicative)
    @parameters begin
        # clouds
        (α_a = 0.275),  [description = "albedo of atmosphere without clouds"]
        (α_s = 0.05), [description = "albedo of surface"]
    end
    if version == :multiplicative
        return α ~ 1 - (1 - α_s)*(1 - α_a)*(1 - α_C)
    elseif version == :additive
        return α ~ α_s + α_a + α_C
    elseif version isa Number
        return α ~ version
    else
        error("incorrect albedo version")
    end
end