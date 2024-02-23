# note; all variables that do not have a process here
# become parameters via `ParameterProcess` by default.
DEFAULT_CCM_PROCESSES = [
    BasicRadiationBalance(),
    # shortwave
    ASR ~ S*(1 - α)*solar_constant,
    IceAlbedoFeedback(),
    DirectAlbedoAddition(),  # Albedo uses fact that cloud albedo is defined as additive
    S ~ 1, # don't make insolation a parameter by default
    # longwave
    LinearOLR(),
    CO2Forcing(), # note that for default CO2 values this is zero forcing
    ParameterProcess(CO2),
    AbsoluteHumidityIGLCRH(),
]
