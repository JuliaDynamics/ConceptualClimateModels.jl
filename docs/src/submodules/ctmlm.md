# [CloudToppedMixedLayerModel](@id ctmlm_docs)

```@docs
CloudToppedMixedLayerModel
```

## Mixed layer

```@docs
CloudToppedMixedLayerModel.mlm_dynamic
CloudToppedMixedLayerModel.bbl_stevens2006_steadystate
CloudToppedMixedLayerModel.temperature_exact
CloudToppedMixedLayerModel.entrainment_velocity
CloudToppedMixedLayerModel.q_liquid
CloudToppedMixedLayerModel.q_saturation
CloudToppedMixedLayerModel.pressure
CloudToppedMixedLayerModel.sst_dynamic
```

## Clouds and decoupling

```@docs
CloudToppedMixedLayerModel.cf_dynamic
CloudToppedMixedLayerModel.cloud_emission_temperature
CloudToppedMixedLayerModel.cloud_emissivity
CloudToppedMixedLayerModel.cloud_albedo
CloudToppedMixedLayerModel.cloud_base_height
CloudToppedMixedLayerModel.decoupling_variable
CloudToppedMixedLayerModel.decoupling_ratios
```

## Radiation

```@docs
CloudToppedMixedLayerModel.cloud_shortwave_warming
CloudToppedMixedLayerModel.cloud_longwave_cooling
CloudToppedMixedLayerModel.mlm_radiative_cooling
CloudToppedMixedLayerModel.downwards_longwave_radiation
CloudToppedMixedLayerModel.albedo
CloudToppedMixedLayerModel.matsunobu_emissivity
```

## Free troposphere

```@docs
CloudToppedMixedLayerModel.free_troposphere_emission_temperature
CloudToppedMixedLayerModel.mlm_s₊
CloudToppedMixedLayerModel.mlm_q₊
```

## Default processes

```@example MAIN
using ConceptualClimateModels
import ConceptualClimateModels.CloudToppedMixedLayerModel as CTMLM
default_processes_eqs(CTMLM)
```
