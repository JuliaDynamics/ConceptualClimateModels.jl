# ConceptualClimateModels.jl

```@docs
ConceptualClimateModels
```

## [Premade symbolic variable instances](@id global_vars)

For convenience, ConceptualClimateModels.jl defines and exports some symbolic variables
that we [list below](@ref list_vars). These are used throughout the library as the
default variables in [implemented processes](@id processes).
When going through documentation strings of processes, such as [`BasicRadiationBalance`](@ref),
you will notice that the function call signature is like:

```julia
BasicRadiationBalance(; T=T, kwargs...)
```

This `T=T` means that the keyword argument `T`, which represents the
"temperature variable", takes the value of `ConceptualClimateModels.T`,
which itself is a premade symbolic variable instance that is exported by
ConceptualClimateModels.jl. You can pass in your own variables instead, by doing
```julia
@variables begin
    (T1_tropics(t) = 290.0), [bounds = (200.0, 350.0), description = "temperature in tropical box 1, in Kelvin"]
end
BasicRadiationBalance(; T=T1_tropics, kwargs...)
```
_(you would also need to give `T1_tropics` to all other processes that utilize temperature)_

Defining variables with the extra `bounds, description` annotations is
useful for integrating with the rest of the functionality of the library.

### [List of premade variables](@id list_vars)

```@example MAIN
using ConceptualClimateModels
PREDEFINED_CCM_VARIABLES
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

# ## Utilities

```@docs
TanhProcess
```

