"""
    GlobalMeanEBM

Submodule of `ConceptualClimateModels` which provides processes
useful in creating global mean energy balance models.
It is inspired by Budyko-Sellers-Ghil type models (without the latitudinal dependence).
The mean does not have to be global, it can be hemispheric, or for other smaller regions.
"""
module GlobalMeanEBM
using ConceptualClimateModels
import NaNMath

include("variables.jl")

# include all processes
for (root, dirs, files) in walkdir(joinpath(@__DIR__))
    for file in files
        file in ("variables.jl", "GlobalMeanEBM.jl") && continue
        include(joinpath(root, file))
    end
end

function __init__()
    # note; all variables that do not have a process here
    # become parameters via `ParameterProcess` by default.
    register_default_process!.(
        [
            BasicRadiationBalance(),
            # shortwave
            ASR ~ S*(1 - α)*solar_constant,
            IceAlbedoFeedback(),
            DirectAlbedoAddition(), # Albedo uses fact that cloud albedo is defined as additive
            # longwave
            BudykoOLR(),
            CO2Forcing(), # for default CO2 values this is zero forcing
            AbsoluteHumidityIGLCRH(),
            # misc
            ΔTLinearRelaxation(),
        ],
        Ref(GlobalMeanEBM)
    )
end


end # Module

export GlobalMeanEBM