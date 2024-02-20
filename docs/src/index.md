# EnergyBalanceModels.jl

## [Premade symbolic variable instances](@id global_vars)

For convenience, EnergyBalanceModels.jl defines and exports some symbolic variables. These are used throughout the library as the default variables in [implemented processes](@id processes).
When going through documentation strings of processes, such as [`BasicRadiationBalance`](@ref),
you will notice that the function call signature is like:

```julia
BasicRadiationBalance(; T=T, kwargs...)
```

This `T=T` means that the keyword argument `T`, which represents the "temperature variable", takes the value of `EnergyBalanceModels.T`, which itself is a premade symbolic variable instance that is exported by EnergyBalanceModels.jl. You can pass in your own variables instead, by doing
```julia
@variables T1_tropics = 310.0
BasicRadiationBalance(; T=T1_tropics, kwargs...)
```

The exported symbolic variables are:

```@example MAIN
using EnergyBalanceModels
PREDEFINED_EBM_VARIABLES
```

## Default values, limits, etc.

All [exported symbolic variable instances](@ref) are defiled with a default value and have plausible physical limits that can be obtained with the following function:

```@docs
physically_plausible_limits(::String)
```

e.g.,

```@example MAIN
physically_plausible_limits(T)
```


## Integration with DynamicalSystems.jl

EnergyBalanceModels.jl integrates with [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/dev/) in many ways.


```@docs
physically_plausible_limits(::DynamicalSystem)
physically_plausible_ic_sampler
physically_plausible_grid
```

# ## Utilities

```@docs
TanhProcess
```

