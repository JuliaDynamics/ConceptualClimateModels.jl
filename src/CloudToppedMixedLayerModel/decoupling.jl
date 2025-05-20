"""
    decoupling_variable(version = :Bretherton1997)

Provide an equation for ``\\mathcal{D}``, the decoupling variable.
"""
function decoupling_variable(version = :Bretherton1997)
    if version == :Bretherton1997
        return 𝒟 ~ RCT*LHF/CTRC
    elseif dversion == :Chung2012
        return 𝒟 ~ LHF/CTRC
    else
        error("unknown")
    end
end

"""
    decoupling_ratios()

Return equations for ``\\mathcal{d}_q, \\mathcal{d}_s`` as in [Datseris2025](@cite).
"""
function decoupling_ratios()
    # Clamping this makes the simulation more stable;
    # the data from de Roode 2016 never exceed 0.5
    α ~ clamp((CTBBL.z_b*CTBBL.RCT/2750)^1.3, 0, 0.5),
    return [𝒹_q ~ α, 𝒹_s ~ 0.5α]
end

# This is a nice plot that creates the same approximate data as in Fig 6 of de Roode 2016:
# fig, axs = axesgrid(2,1; sharex = true, xlabels = "z_cl", ylabels = ["𝒹_q", "𝒹_s"])
# zs = 0:0.1:1.5 # in Km
# 𝒹_q = @. (zs/2.75)^1.3
# lines!(axs[1], zs, 𝒹_q)
# 𝒹_s = @. 0.5*(zs/2.75)^1.3
# lines!(axs[2], zs, 𝒹_s)
# lines!(axs[2], zs, 𝒹_q)
# fig
