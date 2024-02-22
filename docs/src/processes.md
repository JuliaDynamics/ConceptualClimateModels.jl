# [Predefined processes](@id predefined_processes)

Predefined processes are provided in this page.
Those extracted by the literature cite their according resource.
Note that by default all processes utilize the globally-exported [predefined variables](@ref global_vars) of ConceptualClimateModels.jl.

## Temperature

```@docs
BasicRadiationBalance
```

## Temperature difference

```@docs
ΔTLinearRelaxation
ΔTStommelModel
```

## Longwave radiation

```@docs
LinearOLR
LinearClearSkyOLR
EmissivityStefanBoltzmanOLR
EmissivityFeedbackTanh
EmissivitySellers1969
SoedergrenClearSkyEmissivity
```

## Shortwave radiation

```@docs
DirectAlbedoAddition
CoAlbedoProduct
SeparatedClearAllSkyAlbedo
```

## Ice/snow

```@docs
IceAlbedoFeedback
```

## Insolation

```@docs
AstronomicalForcingDeSaedeleer
```

## Forcings

```@docs
CO2Forcing
```

## Clouds

```@docs
CloudAlbedoExponential
CloudAlbedoLinear
BudykoOLR
```

## [Generic processes](@id generic_processes)

Processes that do not depend on any particular physical concept and instead provide
a simple way to create new processes for a given climate variable:

```@docs
ParameterProcess
TimeDerivative
ExpRelaxation
TanhProcess
```

## Default processes

The list of default processes that are used by default in [`processes_to_coupledodes`](@ref) if one does not explicitly provide a list of default processes are:

```@example MAIN
DEFAULT_CCM_PROCESSES
```

