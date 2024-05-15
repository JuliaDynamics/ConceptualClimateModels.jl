module GlobalMeanEBM
using ConceptualClimateModels

include("variables.jl")

# include all processes first
for (root, dirs, files) in walkdir(joinpath(@__DIR__))
    for file in files
        file == "GlobalMeanEBM.jl" && continue
        file == "variables.jl" && continue
        include(joinpath(root, file))
    end
end



end # Module


# Default processes: we need to do this outsime module scope.
# I don't really understand why, but if this code is within the
# module, the processes are not actually registered. :S

# note; all variables that do not have a process here
# become parameters via `ParameterProcess` by default.
let
    using .GlobalMeanEBM
    register_default_process!.(
        [
            BasicRadiationBalance(),
            # shortwave
            ASR ~ S*(1 - α)*GlobalMeanEBM.solar_constant,
            IceAlbedoFeedback(),
            DirectAlbedoAddition(), # Albedo uses fact that cloud albedo is defined as additive
            S ~ 1, # don't make insolation a parameter by default
            # longwave
            BudykoOLR(),
            CO2Forcing(), # for default CO2 values this is zero forcing
            ParameterProcess(CO2),
            AbsoluteHumidityIGLCRH(),
            # misc
            C ~ default_value(C),
            ΔTLinearRelaxation(),
        ],
        Ref(GlobalMeanEBM)
        )

end


