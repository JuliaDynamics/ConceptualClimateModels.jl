module GlobalMeanEBM
using ConceptualClimateModels

include("variables.jl")

# include all processes first
for (root, dirs, files) in walkdir(joinpath(@__DIR__))
    for file in files
        file in ("variables.jl", "GlobalMeanEBM.jl") && continue
        include(joinpath(root, file))
    end
end


# Default processes: we need to do this outsime module scope.
# I don't really understand why, but if this code is within the
# module, the processes are not actually registered. :S

# note; all variables that do not have a process here
# become parameters via `ParameterProcess` by default.

# import .GlobalMeanEBM

function __init__()

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