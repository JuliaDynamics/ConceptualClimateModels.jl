# include all processes first
for (root, dirs, files) in walkdir(joinpath(@__DIR__, "processes"))
    for file in files
        include(joinpath(root, file))
    end
end

# This depends on the `variables.jl` file!
DEFAULT_PROCESSES = [
    BasicRadiationBalance(),
    ASR ~ S*(1 - α)*solar_constant,
    LinearOLR(),
    ParameterProcess(CO2),
    f ~ 0, # could be CO2Forcing() instead
    S ~ 1, # don't make insolation a parameter by default
    # Albedo
    DirectAlbedoAddition(),
    # Default C and BudykoOLR() give same OLR balance as current CERES data
    C ~ default_value(C), # don't make clouds parameter by default
]

export DEFAULT_PROCESSES
# ProcessBasedModelling.default_processes() = DEFAULT_PROCESSES