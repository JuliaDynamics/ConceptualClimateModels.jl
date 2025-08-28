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
            T_lcl ~ s_b - g*CLT/cₚ, # analytically resolve coz no liquid water
            q_x ~ 0,
            s_x ~ 0,
            bbl_emission_temperature(),
            ΔF_q ~ 0.0,

            # Exchanges
            ParameterProcess(SST_X, 0),
            ParameterProcess(S, 400.18),
            q₀ ~ q_saturation(SST), # surface = 100% humid since its ocean
            s₀ ~ SST,
            RH_b ~ q_b/q₀,
            V ~ U*d_c,
            LHF ~ -ρ₀*(ℓ_v/1e3)*V*Δ₀q, # defined as positive, and q is in g/kg
            SHF ~ -ρ₀*V*Δ₀s*cₚ, # defined as positive and s is in cp units

            # Radiation
            L_c ~ σ_SB*ε_C*T_C^4,
            L_b ~ σ_SB*ε_b*T_b^4,
            L_FTR ~ σ_SB*ε_FTR*T_FTR^4,
            Lnet ~ L₀ - Ld, # Ld is whatever reaches the surface
            ASW ~ S*(1 - α) - CRCsw, # Absorbed shortwave at the surface
            ε_FTR ~ matsunobu_emissivity(RH₊, T_FTR),
            ε_b ~ matsunobu_emissivity(RH_b, T_b),
            L₀ ~ σ_SB*SST^4,
            albedo(),

            # Clouds
            z_cb ~ z_lcl,
            z_ct ~ z_b,
            CLT ~ max(z_ct - z_cb, 0),
            RCT ~ min(CLT/z_b, 1), # note: not devided by cloud top but by boundary layer.
            CRC ~ CRClw - CRCsw,
            CTRC ~ CTRClw - CRCsw,
            LWP ~ liquid_water_path_linear(),
            cloud_emission_temperature(),

            # Decoupling
            λ_s ~ 0,
            λ_q ~ 0,
        ],
        Ref(CloudToppedMixedLayerModel)
    )
end
