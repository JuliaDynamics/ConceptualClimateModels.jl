"""
    free_troposphere_emission_temperature(γ = 1.0; add_co2 = true)

Return an equation for ``T_{FTR}``. `γ` is as in [Datseris2025](@cite).
`add_co2` will add an additional warming term `ECS_CO2*log2(CO2/400)`.
"""
function free_troposphere_emission_temperature(version = 1.0; add_co2 = true)
    # We measure the subtropical component of T_FTR from T₊
    # not from T_t, because that's what makes sense due to physical continuity.
    # T₊ already has a CO2 component so we are careful here to add it only
    # to the tropical component of T_FTR!
    @parameters begin
        (T_FTR_0 = 290.0), [description = "prescribed temperature emission of free troposphere without CO2 effects, K"]
        (T₊_0 = 5.0), [description = "Subtracted temperature from T₊ when defining T_FTR as a T₊, K"]
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
    eq = weight*T_FTR_0 + (1 - weight)*(T₊ - T₊_0) - δ_FTR
    if add_co2 == true
        eq += ECS_CO2*log2(CO2/400)
    end
    return T_FTR ~ eq
end

###########################################################################################
# coupling BBL with free troposphere
###########################################################################################
"""
    bbl_s₊(
        version = :difference;
        cloud_effect = false,
        CO2_effect = false,
    )

Provide equation for ``s_+``. To do this, a boundary condition must be provided
that is a fixed parameter. `version` argument decides this:

- `:difference`: the starting temperature difference across inversion is a fixed parameter.
- `:temperature`: the starting temperature after the inversion is a fixed parameter.
- `:static_energy`: the starting moist static energy after the inversion is a fixed parameter.

Besides these, we can also specify whether CO2 increase also increases temperature difference,
and whether decreasing ``C`` decreases temperature difference due to cloud thinning
as in [Singer2023a](@cite).
"""
function bbl_s₊(
        inversion_fixing = :difference;
        cloud_effect = false,
        CO2_effect = false,
    )
    @parameters begin
        (Δ₊T_C = 10.0), [description = "temperature decrease in the inversion due to cloud thinning (max lost K for 100% C.F.), K"]
        (s₊_0 = 300.0), [description = "prescribed moist static energy above inversion normalized by cₚ, K"]
        (T₊_0 = 292.0), [description = "prescribed temperature above inversion without CO2 or cloud effects, K"]
    end

    # First, prepare the augmentation of the inversion
    Δ₊T_aux = δ_Δ₊T
    if cloud_effect
        # This is done in Singer & Schendier but I don't think it is correct.
        # It is not justified in any way that I found reasonable.
        Δ₊T_aux += - Δ₊T_C*(1 - C)
    end
    if CO2_effect
        Δ₊T_aux += ECS_CO2*log2(CO2/400)
    end

    # Now we see which type of inversion we keep prescribed
    if inversion_fixing == :difference
        eqs = [
            s₊ ~ T₊ + g*z_b/cₚ,
            T₊ ~ T_t + Δ₊T,
            Δ₊T ~ Δ₊T_aux
        ]
    elseif inversion_fixing == :temperature
        eqs = [
            T₊ ~ T₊_0 + Δ₊T_aux,
            s₊ ~ T₊ + g*z_b/cₚ,
            Δ₊T ~ T₊ - T_t, # observable
        ]
    elseif fixed == :static_energy
        eqs = [
            s₊ ~ s₊_0 + Δ₊T_aux,
            T₊ ~ s₊ - g*z_b/cₚ,
            Δ₊T ~ T₊ - T_t, # observable
        ]
    else
        error("incorrect specification for what to stay fixed!")
    end
    return eqs
end

"""
    bbl_q₊(version = :relative)

Provide equation for ``q_+``. If `version = :relative` then
make free tropospheric relative humidity a free parameter.
Else if `version = :constant` then make ``q_+`` itself a parameter.
"""
function bbl_q₊(humidity_fixing = :relative)
    if humidity_fixing == :relative
        return q₊ ~ RH₊ * q_saturation(T₊)
    elseif humidity_fixing == :constant
        return ParameterProcess(q₊, 1.5)
    end
end

function bbl_Δ₊sᵥ() # q_b are global symbols
    # definition of s:
    # s = T + g*z/cₚ - (ℓ_v/1e3)*q_liquid(T, q)/cₚ # units of g/kg for `q` hence /1e3
    # divide by c_p, and no liquid after inversion
    ql = q_liquid(T_t, q_b, z_b) # liquid water at top of cloud
    Δ₊Tᵥ = T₊*(1 + 0.608q₊/1e3) - T_t*(1 + 0.608q_b/1e3 - 1.608ql/1e3)
    return Δ₊sᵥ ~ Δ₊Tᵥ - (ℓ_v/1e3)*ql/cₚ
end
