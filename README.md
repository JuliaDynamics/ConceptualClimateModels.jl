# ConceptualClimateModels.jl

[![docsdev](https://img.shields.io/badge/docs-dev-lightblue.svg)](https://juliadynamics.github.io/ProcessBasedModelling,jl/dev/)
[![docsstable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliadynamics.github.io/ConceptualClimateModels.jl/stable/)
[![CI](https://github.com/JuliaDynamics/ConceptualClimateModels.jl/workflows/CI/badge.svg)](https://github.com/JuliaDynamics/ConceptualClimateModels.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/JuliaDynamics/ConceptualClimateModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaDynamics/ConceptualClimateModels.jl)
[![Package Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/ProcessBasedModelling)](https://pkgs.genieframework.com?packages=ProcessBasedModelling)

ConceptualClimateModels.jl is a Julia package for creating and analysing conceptual
models of climate, such as energy balance models, glaciation cycle models, or climate tipping models.
Such conceptual models are simplified representation of basic climate components,
and the processes that connect them, such as flows of energy or mass.
Within this context such models are typically coupled ordinary differential
equations (with partial or stochastic DEs also being possible).

ConceptualClimateModels.jl accelerates both modelling and analysis aspects of working
with such models by:

- Building upon [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/)
  for creating equations from _symbolic expressions_.
- Utilizing [ProcessBasedModelling.jl](https://github.com/JuliaDynamics/ProcessBasedModelling.jl?tab=readme-ov-file)
  to provide a field-specific framework that allows easily testing different physical
  hypotheses regarding how climate variables couple to each
  other, or how climate processes are represented by equations.
- Offering many predefined processes from current literature and ongoing research.
- Being easy to extend with more climate variables or physical processes.
- Allowing the straightforward coupling of different conceptual models with each other.
- Automating the naming of custom parameters relating to existing climate processes.
- Integrating with [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/dev/)
  to automate the start-up phase of typical nonlinear dynamics based workflows.

with other features coming soon, such as:

- Support for latitudinal models (where each variable is vector-valued over latitude circles)
- Support for stochastic ordinary differential equations

To install it, run `import Pkg; Pkg.add("ConceptualClimateModels")`.

All further information is provided in the documentation, which you can either find
[online](https://juliadynamics.github.io/ConceptualClimateModels.jl/stable/) or build
locally by running the `docs/make.jl` file.

ConceptualClimateModels.jl development is funded by UKRI's Engineering and Physical Sciences Research Council, grant no. EP/Y01653X/1 (grant agreement for a EU Marie Sklodowska-Curie Postdoctoral Fellowship).