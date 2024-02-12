# Predefined processes

Predefined processes are provided in this page.
Those extracted by the literature cite their according resource.
Note that by default all processes utilize the globally-exported [predefined variables](@ref variables) of EnergyBalanceModels.jl.
This is conveyed in the documentation strings by writing something like
```
f(; T = T, kw...)
```
this means that the variable representing "T" will be by default the global symbolic variable `T`. You could of course provide any other variable, such as `T1` or `T2` if you have two "boxes" with different temperatures for example.

## Temperature

```@docs
BasicRadiationBalance
ΔTLinearRelaxation
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
EmissivityCumulativeSubtraction
BudykoOLR
```
