# Examples

In this page we go through various examples of combining processes to make
models that have been already used in the literature, or using DynamicalSystems.jl
or other packages to analyse conceptual climate models.

## Insolation and Snowball Earth Bistability

The origins of energy balance models ([Sellers1969, Budyko1969, Ghil1981, North1981](@cite))
examined the impact of variations in insolation on the global climate.
In particular, they studied how a global-mean temperature energy balance model
yielded bi-stable hysterisis between a cold "snowball" state and a warm Earth,
as the solar constant was varied.

We can easily replicate this behaviour and create the same model as [Budyko1969](@cite),
by combining the processes:

```@example MAIN
budyko_processes = [
    BasicRadiationBalance(),
    IceAlbedoFeedback(),
    BudykoOLR(),
    C ~ default_value(C), # make cloud fraction a fixed constant
]

budyko = model
```