# GlobalMeanEBM

```@docs
GlobalMeanEBM
```

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

## Water vapor

```@docs
saturation_vapor_pressure
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
AdditionProcess
TanhProcess
```

## Default processes

```@example MAIN
using ConceptualClimateModels.GlobalMeanEBM
default_processes_eqs(GlobalMeanEBM)
```
