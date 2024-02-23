# Examples

In this page we go through various examples of combining processes to make
models that have been already used in the literature, or using DynamicalSystems.jl
or other packages to analyse conceptual climate models.

## Insolation and Snowball Earth Bistability

The origins of energy balance models ([Sellers1969, Budyko1969, Ghil1981, North1981](@cite))
examined the impact of variations in insolation on the global climate.
In particular, they studied how simple energy balance models with only ice-albedo and water vapor feedbacks yielded bi-stable hysteresis between a cold "snowball" state and
a warm Earth, as the solar constant was varied.

We can easily replicate this behaviour by creating a
global-mean temperature model, by combining the processes:

```@example MAIN
using ConceptualClimateModels

budyko_processes = [
    BasicRadiationBalance(),
    LinearOLR(),
    IceAlbedoFeedback(; min = 0.3, max = 0.7),
    α ~ α_ice,
    ParameterProcess(S), # insolation is a parameter
    f ~ 0, # no external forcing
    # absorbed solar radiation has a default process
]

budyko = processes_to_coupledodes(budyko_processes)
dynamical_system_summary(budyko)
```

The water vapor feedback is encapsulated in the form of the OLR, which is linear
due to water vapor [Koll2018](@cite).
