"Return τSST * d(SST)/dt = ASW - Lnet - LHF - SHF + SST_X ."
function sst_dynamic()
    TimeDerivative(SST, (ASW - Lnet - LHF - SHF + SST_X), τ_SST)
end

function __init__()
    register_default_process!.(
        Union{Process, Equation}[
            # Boundary layer
            ρ₀ ~ p₀/(Rd*SST), # moist air surface density in kg/m³
            Δ₊s ~ s₊ - s_b,
            Δ₊q ~ q₊ - q_b,
            Δ₀q ~ q_b - q₀,
            Δ₀s ~ s_b - s₀,
            ζ ~ 0,
            w_m ~ 0,
            w_v ~ 0,
            T_t ~ temperature_exact(z_b, s_b, q_b), # it isn't guaranteed that this will be used!
            T_lcl ~ s_b - g*CLT*z_b/cₚ, # from the cloud base height exact function
            q_x ~ 0,
            s_x ~ 0,
            bbl_emission_temperature(),

            # Exchanges
            ParameterProcess(SST_X, 0),
            ParameterProcess(S, 400.18),
            q₀ ~ q_saturation(SST), # surface = 100% humid since its ocean
            s₀ ~ SST,
            RH_b ~ q_b/q₀,
            V ~ U*c_d,
            LHF ~ -ρ₀*(ℓ_v/1e3)*V*Δ₀q, # defined as positive, and q is in g/kg
            SHF ~ -ρ₀*V*Δ₀s*cₚ, # defined as positive and s is in cp units

            # Radiation
            L_c ~ σ_SB*ε_c*T_c^4,
            L_b ~ σ_SB*ε_b*T_b^4,
            L_FTR ~ σ_SB*ε_FTR*T_FTR^4,
            Lnet ~ L₀ - Ld, # Ld is whatever reaches the surface
            ASW ~ S*(1 - α) - CRCsw, # Absorbed shortwave at the surface
            ε_FTR ~ matsunobu_emissivity(RH₊, T_FTR),
            ε_b ~ matsunobu_emissivity(RH_b, T_b),
            L₀ ~ σ_SB*SST^4,
            albedo(),

            # Clouds
            CRC ~ CRClw - CRCsw,
            CTRC ~ CTRClw - CRCsw,
            LWP ~ liquid_water_path(T_t, CLT, z_b, s_b, q_b),
            cloud_emission_temperature(),

            # Decoupling
            𝒹_s ~ 0,
            𝒹_q ~ 0,
        ],
        Ref(CloudToppedMixedLayerModel)
    )
end
