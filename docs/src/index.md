# ConceptualClimateModels.jl

```@docs
ConceptualClimateModels
```


## Default values, limits, etc.

All [exported symbolic variable instances](@ref) are defiled with a default value and have plausible physical limits that can be obtained with the following function:

```@docs
physically_plausible_limits(::Any)
```

e.g.,

```@example MAIN
physically_plausible_limits(T)
```


## Integration with DynamicalSystems.jl

ConceptualClimateModels.jl integrates with [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/dev/) in many ways.


```@docs
physically_plausible_limits(::DynamicalSystem)
physically_plausible_ic_sampler
physically_plausible_grid
```

