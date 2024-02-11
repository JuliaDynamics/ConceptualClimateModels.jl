# EnergyBalanceModels.jl

Any quantity that may potentially be a dynamic (state) variable, or an explicit function of time, is defined in `src/variables.jl` and must have a default value when defined. They should also be exported and in the rest of the source code they are accessed from the module-level scope. _TODO: Example here_.
Variables that may be dynamic (state) variables must also obtain limits in the `physicall_plausible_limits` function which is in the same file.
Variables that can never be dynamic but are guaranteed to be observed, such as the outgoing longwave radiation, should be preferably defined in the same file but without a default value or limits.

# TODO:
explain that all processes by default use the global variables,
which is why docstrings are written this way.

EnergyBalanceModels

"""
module defining various equations and models representing energy balance
processes. It utilizes ModelingToolkit.jl to create a flexible framework
for adding or remove dynamical state variables and adding or removing
equations that represent various processes in the earth system.
"""