# Predefined processes

Predefined processes are provided in this page.
Those extracted by the literature cite their according resource.
Note that by default all processes utilize the globally-exported [predefined variables](@ref global_vars) of ConceptualClimateModels.jl.

## Temperature

```@docs
BasicRadiationBalance
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
