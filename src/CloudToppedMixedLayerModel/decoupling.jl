"""
    decoupling_variable(version = :Bretherton1997)

Provide an equation for ``\\Lambda``, the decoupling variable.
"""
function decoupling_variable(version = :Bretherton1997)
    if version == :Bretherton1997
        return Λ ~ RCT*LHF/CTRC
    elseif dversion == :Chung2012
        return Λ ~ LHF/CTRC
    else
        error("unknown")
    end
end

"""
    decoupling_ratios()

Return equations for ``\\lambda_q, \\lambda_s`` as in [Datseris2025](@cite).
"""
function decoupling_ratios()
    # Clamping this makes the simulation more stable;
    # the data from de Roode 2016 never exceed 0.5
    α ~ clamp((CTBBL.z_b*CTBBL.RCT/2750)^1.3, 0, 0.5),
    return [λ_q ~ α, λ_s ~ 0.5α]
end

# This is a nice plot that creates the same approximate data as in Fig 6 of de Roode 2016:
# fig, axs = axesgrid(2,1; sharex = true, xlabels = "z_cl", ylabels = ["λ_q", "λ_s"])
# zs = 0:0.1:1.5 # in Km
# λ_q = @. (zs/2.75)^1.3
# lines!(axs[1], zs, λ_q)
# λ_s = @. 0.5*(zs/2.75)^1.3
# lines!(axs[2], zs, λ_s)
# lines!(axs[2], zs, λ_q)
# fig
