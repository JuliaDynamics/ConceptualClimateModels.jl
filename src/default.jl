# include all processes first
for (root, dirs, files) in walkdir(joinpath(@__DIR__, "processes"))
    for file in files
        include(joinpath(root, file))
    end
end

# This depends on the `variables.jl` file!
DEFAULT_CCM_PROCESSES = [
    BasicRadiationBalance(),
    ASR ~ S*(1 - α)*solar_constant,
    LinearOLR(),
    CO2Forcing(), # note that for default CO2 values this is zero forcing
    ParameterProcess(CO2),
    S ~ 1, # don't make insolation a parameter by default
    # Albedo
    DirectAlbedoAddition(),
]

export DEFAULT_CCM_PROCESSES