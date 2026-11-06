"""
    free_troposphere_emission_temperature(γ = 1.0; add_co2 = true)

Return an equation for ``T_{FTR}``. `γ` is as in [Datseris2026](@cite).
`add_co2` will add an additional warming term `ECS_CO2*log2(CO2/400)`.
Introduces two extra parameters, `T_FTR_0, Tᵤ_0`.
"""
function free_troposphere_emission_temperature(version = 1.0; add_co2 = true, Tᵤ_0_init = 5.0)
    # We measure the subtropical component of T_FTR from Tᵤ
    # not from T_t, because that's what makes sense due to physical continuity.
    # Tᵤ already has a CO2 component so we are careful here to add it only
    # to the tropical component of T_FTR!
    @parameters begin
        (T_FTR_0 = 290.0), [description = "prescribed temperature emission of free troposphere without CO2 effects, K"]
        (Tᵤ_0 = Tᵤ_0_init), [description = "Subtracted temperature from Tᵤ when defining T_FTR as a Tᵤ, K"]
    end
    if version == :strong # strong influence of subtropical
        weight = 0.0
    elseif version == :weak # some influence
        weight = 0.5
    elseif version == :none # only tropical
        weight = 1.0
    elseif version isa Number
        weight = version
    else
        error("incorret gradient for T_FTR")
    end
    eq = weight*T_FTR_0 + (1 - weight)*(Tᵤ - Tᵤ_0) - δ_FTR
    if add_co2 == true
        eq += ECS_CO2*log2(CO2/400)
    end
    return T_FTR ~ eq
end

###########################################################################################
# coupling BBL with free troposphere
###########################################################################################
"""
    mlm_sᵤ(
        version = :difference;
        cloud_effect = false,
        CO2_effect = false,
    )

Provide equation for ``s_+`` depending on `version` with options:

- `:difference`: the temperature difference across inversion is a fixed parameter.
- `:temperature`: the temperature after the inversion is a fixed parameter.
- `:static_energy`: the moist static energy after the inversion is a fixed parameter.
- `:lapse_rate`: the temperature after the inversion is fixed and given by a prescribed
  lapse rate, ``T_+ = T_{+,ref} + \\Gamma_T (z_b - 1000)`` as in [Salazar2023](@cite)
  which introduces two additional parameters: `Γ_T = 6.5e-3, Tᵤ_ref = 290.0`

Besides these, we can also specify whether CO2 increase also increases temperature difference,
and whether decreasing ``C`` decreases temperature difference due to cloud thinning
as in [Singer2023a](@cite).
"""
function mlm_sᵤ(
        inversion_fixing = :difference;
        cloud_effect = false,
        CO2_effect = false,
    )
    @parameters begin
        (ΔᵤT_C = 10.0), [description = "temperature decrease in the inversion due to cloud thinning (max lost K for 100% C.F.), K"]
        (sᵤ_0 = 300.0), [description = "prescribed moist static energy above inversion normalized by cₚ, K"]
        (Tᵤ_0 = 292.0), [description = "prescribed temperature above inversion without CO2 or cloud effects, K"]
        (Tᵤ_ref = 290.0), [description = "reference temperature above inversion (for lapse rate), K"]
        (Γ_T = 6.5e-3), [description = "free tropospheric temperature lapse rate, K/m"]
    end

    # First, prepare the augmentation of the inversion
    ΔᵤT_aux = δ_ΔᵤT
    if cloud_effect
        # This is done in Singer & Schneider but I don't think it is correct.
        # It is not justified in any way that I found reasonable.
        ΔᵤT_aux += - ΔᵤT_C*(1 - C)
    end
    if CO2_effect
        ΔᵤT_aux += ECS_CO2*log2(CO2/400)
    end

    # Now we see which type of inversion we keep prescribed
    if inversion_fixing == :difference
        eqs = [
            sᵤ ~ Tᵤ + g*z_b/cₚ,
            Tᵤ ~ T_t + ΔᵤT,
            ΔᵤT ~ ΔᵤT_aux
        ]
    elseif inversion_fixing == :temperature
        eqs = [
            Tᵤ ~ Tᵤ_0 + ΔᵤT_aux,
            sᵤ ~ Tᵤ + g*z_b/cₚ,
        ]
    elseif inversion_fixing == :static_energy
        eqs = [
            sᵤ ~ sᵤ_0 + ΔᵤT_aux,
            Tᵤ ~ sᵤ - g*z_b/cₚ,
        ]
    elseif inversion_fixing == :lapse_rate
        eqs = [
            Tᵤ ~ Tᵤ_ref - Γ_T*(z_b - 1000.0),
            sᵤ ~ Tᵤ + g*z_b/cₚ,
        ]
    else
        error("incorrect specification for what to stay fixed!")
    end
    return eqs
end

"""
    mlm_qᵤ(version = :relative)

Provide equation for ``q_+``. If `version = :relative` then
make free tropospheric relative humidity `RHᵤ` a free parameter.
Else if `version = :constant` then make `qᵤ` itself a parameter.
Else if `version = :lapse_rate` prescribe ``q_+ = q_{+, ref} - \\Gamma_q(z_b - 1000)``
which introduces parameters `Γ_q = 1.5e-3, qᵤ_ref = 2.0`, inspired by [Park2004](@cite).
"""
function mlm_qᵤ(humidity_fixing = :relative)
    if humidity_fixing == :relative
        return qᵤ ~ RHᵤ * q_saturation(Tᵤ)
    elseif humidity_fixing == :constant
        return ParameterProcess(qᵤ, 1.5)
    elseif humidity_fixing == :lapse_rate
        @parameters begin
            (qᵤ_ref = 2.0), [description = "reference specific humidity above inversion (for lapse rate), q/kg"]
            (Γ_q = 1.5e-3), [description = "specific humidity lapse rate, g/kg/m"]
        end
        return qᵤ ~ qᵤ_ref - Γ_q * (z_b - 1000.0)
    end
end

function bbl_Δᵤsᵥ() # q_b are global symbols
    # definition of s:
    # s = T + g*z/cₚ - (ℓ_v/1e3)*q_liquid(T, q)/cₚ # units of g/kg for `q` hence /1e3
    # divide by c_p, and no liquid after inversion
    ql = q_liquid(T_t, q_b, z_b) # liquid water at top of cloud
    ΔᵤTᵥ = Tᵤ*(1 + 0.608qᵤ/1e3) - T_t*(1 + 0.608q_b/1e3 - 1.608ql/1e3)
    return Δᵤsᵥ ~ ΔᵤTᵥ - (ℓ_v/1e3)*ql/cₚ
end
