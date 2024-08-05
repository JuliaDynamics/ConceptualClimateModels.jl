# GlobalMeanEBM

```@docs
GlobalMeanEBM
```

## Default variables

```@example MAIN
using ConceptualClimateModels
GlobalMeanEBM.globalmeanebm_variables
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

## Default processes

```@example MAIN
using ConceptualClimateModels.GlobalMeanEBM
default_processes_eqs(GlobalMeanEBM)
```
