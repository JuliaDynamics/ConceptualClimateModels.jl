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
GlobalMeanEBM.BasicRadiationBalance
```

## Temperature difference

```@docs
GlobalMeanEBM.ΔTLinearRelaxation
GlobalMeanEBM.ΔTStommelModel
```

## Longwave radiation

```@docs
GlobalMeanEBM.LinearOLR
GlobalMeanEBM.LinearClearSkyOLR
GlobalMeanEBM.EmissivityStefanBoltzmanOLR
GlobalMeanEBM.EmissivityFeedbackTanh
GlobalMeanEBM.EmissivitySellers1969
GlobalMeanEBM.SoedergrenClearSkyEmissivity
```

## Shortwave radiation

```@docs
GlobalMeanEBM.DirectAlbedoAddition
GlobalMeanEBM.CoAlbedoProduct
GlobalMeanEBM.SeparatedClearAllSkyAlbedo
```

## Ice/snow

```@docs
GlobalMeanEBM.IceAlbedoFeedback
```

## Water vapor

```@docs
GlobalMeanEBM.saturation_vapor_pressure
```

## Insolation

```@docs
GlobalMeanEBM.AstronomicalForcingDeSaedeleer
```

## Forcings

```@docs
GlobalMeanEBM.CO2Forcing
```

## Clouds

```@docs
GlobalMeanEBM.CloudAlbedoExponential
GlobalMeanEBM.CloudAlbedoLinear
GlobalMeanEBM.BudykoOLR
```

## Default processes

```@example MAIN
import ConceptualClimateModels.GlobalMeanEBM
default_processes_eqs(GlobalMeanEBM)
```
