# Tutorial

With ConceptualClimateModels.jl one makes differential equation systems from _processes_.
A _process_ is simply a particular equation defining a climate variable.
A vector of processes is composed by the user, and given to the main function [`processes_to_coupledodes`](@ref) which bundles them into a system of equations.

!!! note "Familiarity with DynamicalSystems.jl and ModelingToolkit.jl"
    ConceptualClimateModels.jl builds on [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/) for building the equations representing the climate model, and it builds on [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/dev/) to analyze the models. Familiarity with either package is good to have, and will allow you to faster and better understand the concepts discussed here. Nevertheless familiarity is actually optional as the steps required to use ConceptualClimateModels.jl are all explained in this tutorial.


## Introductory example

Let's say that we want to create the most basic energy balance model,
$$
c_T \frac{dT}{dt} = ASR - OLR + f
$$
where $ASR$ is the absorbed solar radiation given by
$$
ASR = S (1-\alpha)
$$
with $\alpha$ the planetary albedo, $OLR$ is the outgoing longwave radiation given by the linearized expression
$$
OLR = A + BT
$$
and $f$ some radiative forcing at the top of the atmosphere, that is based on CO2 concentrations and given by
$$
f = 3.7\log_2\left(\frac{CO_2}{400}\right)
$$
with CO2 concentrations in ppm.

To create this model with ConceptualClimateModels.jl while providing the least information possible we can do:

```@example MAIN
using ConceptualClimateModels

processes = [
    BasicRadiationBalance(),
    LinearOLR(),
    CO2Forcing(), # note that for default CO2 values this is zero forcing
]

ds = processes_to_coupledodes(processes)
```

This is a dynamical system from [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/stable/) that is generated via symbolic expressions based on [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/), utilizing the process-based approach of [ProcessBasedModelling.jl](https://juliadynamics.github.io/ProcessBasedModelling.jl/stable/).

As the dynamical system is made by symbolic expressions, these can be
obtained back at any time:

```@example MAIN
using DynamicalSystems
# access the MTK model that stores the symbolic bindings
mtk = referrenced_sciml_model(ds)
# show the equations of the dynamic state variables of the dynamical system
equations(mtk)
```

```@example MAIN
# show state functions that are observable,
# i.e., they do not have a time derivative, they are not dynamic state variables
observed(mtk)
```

```@example MAIN
# show parameters
parameters(mtk)
```

The symbolic variables can also be used to query or alter the dynamical system,

```@example MAIN
observe_state(ds, T)
```

```@example MAIN
# access symbolic parameter CO2_0 from the tracked symbolic list of the model
set_parameter!(ds, mtk.CO2_0, 800)
```

```@example MAIN
current_parameter(ds, mtk.CO2_0)
```

Let's unpack the steps that led to this level of automation.

## Processes

A process is conceptually an _equation the defines a climate variable or observable_.
All processes that composed the system are then composed into a set of differential equations via [`processes_to_coupledodes`](@ref) (or [`processes_to_mtkmodel`](@ref)) that represent the climate model.

For example, the process
```@example MAIN
T_process = BasicRadiationBalance()
```
is the process defining the variable ``T``, representing temperature. We can learn this by either reading the documentation string of [`BasicRadiationBalance`](@ref), or querying it directly:
```@example MAIN
# This is the equation created by the process
lhs(T_process) ~ rhs(T_process)
```

Notice that this process does not further define e.g. outgoing longwave radiation `OLR(t)`.
That is why in the original example we also provided [`LinearOLR`](@ref), which defines it:
```@example MAIN
OLR_process = LinearOLR()
lhs(OLR_process) ~ rhs(OLR_process)
```

Each physical "observable" or variable that can be configured in the system has its own process.
This allows very easily exchanging the way processes are represented by equations without having to alter many equations. For example, if instead of `LinearOLR` we have provided
```@example MAIN
OLR_process = EmissivityStefanBoltzmanOLR()
lhs(OLR_process) ~ rhs(OLR_process)
```
then we would have used a Stefan-Boltzmann grey-atmosphere representation for the outgoing longwave radiation.

## Default processes

Hold on a minute though, because in the original processes we provided,
```julia
processes = [
    BasicRadiationBalance(),
    LinearOLR(),
    CO2Forcing(),
]
```
there was no process that defined the absorbed solar radiation ``ASR``!

Well, ConceptualClimateModels.jl has a list of [default processes](@ref default_processes) that are automatically included in every call to [`processes_to_coupledodes`](@ref).
The default processes for the ASR is $ASR = S(1-\alpha)$ with $S$ the solar constant.
The function [`processes_to_coupledodes`](@ref) goes through all processes the user provided and identifies variables that themselves do not have a process.
It then checks the list of default processes and attempt to assign one to these variables.

If there are no default processes, it makes the variables themselves parameters with the same name but with a subscript 0.

For example, let's assume that we completely remove default processes like so:

```@example MAIN
processes = [
    BasicRadiationBalance(),
    LinearOLR(),
    CO2Forcing(), # note that for default CO2 values this is zero forcing
    ASR ~ S*(1-α), # add the processes for ASR, but not for S or α
]

# note the empty list as 2nd argument, which is otherwise
# the default processes
mtk = processes_to_mtkmodel(processes, [])

equations(mtk)
```

You will notice that here there is no equation $\alpha = \alpha_bg + \alpha_ice + \alpha_clouds$, which is the default processes assigned to variable $\alpha$.
Instead, the process that was just created for it was $\alpha = \alpha_0$,
where $\alpha_0$ is now a _parameter_ of the system (i.e., it can be altered after creating
the system). The value of $\alpha_0$ is the default value of $\alpha$:
```@example MAIN
default_value(α)
```

```@example MAIN
# current value of the _parameter_ α_0 (not the variable!)
default_value(mtk.α_0)
```

When this automation occurs a warning is thrown:

```
┌ Warning: Variable α was introduced in process of variable ASR(t).
│ However, a process for α was not provided,
│ and there is no default process for it either.
│ Since it has a default value, we make it a parameter by adding a process:
│ `ParameterProcess(α)`.
└ @ ProcessBasedModelling ...\ProcessBasedModelling\src\make.jl:65
```

[`ParameterProcess`](@ref) is the most trivial process: it simply means that the corresponding variable does not have any physical process and rather it is a system
parameter.

This automation does not occur if there is no default value.
For example, variables that can never be dynamic state variables, such as $ASR$, do not have a default value. If we have not assigned a process for $ASR$, the system construction would error instead:

```julia
processes = [
    BasicRadiationBalance(),
    LinearOLR(),
    CO2Forcing(),
]

mtk = processes_to_mtkmodel(processes, [])
```

```
ERROR: ArgumentError: Variable ASR was introduced in process of
variable T(t). However, a process for ASR was not provided,
there is no default process for ASR, and ASR doesn't have a default value.
Please provide a process for variable ASR.
```

These warnings and errors are always perfectly informative.
They tell us exactly which variable does not have a process, and exactly which other process introduced the process-less variable first.
This makes the modelling experience stress-free, especially when large and complex
models are being created.

## Adding more processes

ConceptualClimateModels.jl provides an increasing list of
[predefined processes](@ref predefined_processes) that you can use out of
the box to compose climate models. The predefined processes all come from existing
literature and cite their source via BiBTeX.

It is also very easy to make new processes on your own.
The simplest way to make a process is to just provide an equation for it
with the l.h.s. of the equation being the variable the process defines.

For example:

```@example MAIN
@variables x(t) = 0.5 # all variables must be functions of (t)
x_process = x ~ 0.5*T^2 # x is just a function of temperatur
```

A more re-usable approach however is to create a function that generates a process
or create a new process type as we describe in [making new processes](@ref new_processes).

## [Premade symbolic variable instances](@id global_vars)

You might be wondering, when we wrote the equation `ASR ~ S*(1-α)` for the $ASR$ process,
or when we wrote `x ~ 0.5 * T^2`, where did the variable bindings `ASR, S, α` come from?
For convenience, ConceptualClimateModels.jl defines and exports some symbolic variables
for typical climate quantities that have default processes. We list all of these
[below](@ref list_vars). These default bindings are used throughout the library as the
default variables in [predefined processes](@id predefined_processes).

Crucially, these default bindings are _symbolic variables_. They are defined as
```julia
@variables begin
    T(t) = 0.5 # ...
    # ...
end
```
which means that expressions that involve them result in symbolic expressions,
```@example MAIN
A2 = 0.5
B2 = 0.5
OLR2 = A2 + B2*T
```

In contrast, if we did instead
```@example MAIN
T2 = 0.5 # _not_ symbolic!
OLR2 = A2 + B2*T2
```
This `OLR2` is _not_ a symbolic expression and _cannot_ be used to represent a process.

When going through documentation strings of [predefined processes](@ref predefined_processes),
such as [`BasicRadiationBalance`](@ref),
you will notice that the function call signatures are like:

```julia
BasicRadiationBalance(; T=T, kwargs...)
```

This `T=T` means that the keyword argument `T`, which represents the
"temperature variable", takes the value of `ConceptualClimateModels.T`,
which itself is a premade _symbolic_ variable instance that is exported by
ConceptualClimateModels.jl. You can pass in your own variables instead, by doing
```julia
@variables begin
    (T1_tropics(t) = 290.0), [bounds = (200.0, 350.0), description = "temperature in tropical box 1, in Kelvin"]
end
process = BasicRadiationBalance(; T=T1_tropics, kwargs...)
```
_(you would also need to give `T1_tropics` to all other processes that utilize temperature)_

Defining variables with the extra `bounds, description` annotations is
useful for integrating with the rest of the functionality of the library.

### [List of premade variables](@id list_vars)

All premade variables have a default value, a description, and plausible physical bounds
for integration with DynamicalSystems.jl.

```@example MAIN
using ConceptualClimateModels
PREDEFINED_CCM_VARIABLES
```

## Automatically named parameters

The majority of [predefined processes](@ref predefined_processes) create symbolic
parameters that are automatically named based on the variables governing the processes.
This default behaviour can be altered in two ways.

For example, [`IceAlbedoFeedback`](@ref) adds named parameters to the equations
whose name is derived from the name of the variable representing
ice albedo:
```@example MAIN
@variables my_ice_α(t) = 0.1 # don't forget the `(t)`!
ice_process = IceAlbedoFeedback(; α_ice = my_ice_α)
processes = [ice_process]

mtk = processes_to_mtkmodel(processes)
equations(mtk)
```

```@example MAIN
parameters(mtk)
```

We can alter this behaviour by either providing our own named parameters
to one of the keywords of the process, or wrapping a value around
`LiteralParameter` to replace the parameter by a literal constant, like so:

```@example MAIN
@parameters myfreeze = 260.0
ice_process = IceAlbedoFeedback(;
    α_ice = my_ice_α, Tfreeze = myfreeze, max = LiteralParameter(0.9)
)

mtk = processes_to_mtkmodel([ice_process])
equations(mtk)
```

```@example MAIN
parameters(mtk)
```


## API Reference

```@docs
processes_to_coupledodes
processes_to_mtkmodel
```

### [Making new processes](@id new_processes)

To make a new processes you can:

1. Create a _function_ that given some keyword arguments (including which symbolic
   variables to use) uses one of the existing [generic processes](@ref generic_processes)
   to make and return a process instance. Or, it can return an equation directly,
   provided it satisfies the format of [`processes_to_mtkmodel`](@ref).
   For an example of this, see the source code of [`SeparatedClearAllSkyAlbedo`](@ref)
   or [`EmissivityFeedbackTanh`](@ref).
2. Create a new `Process` subtype. This is preferred, because
   it leads to much better printing/display of the list of processes.
   For an example of this, see the source code of [`IceAlbedoFeedback`](@ref).
   To create a new `Process` see the [API of ProcessBasedModelling.jl](https://juliadynamics.github.io/ProcessBasedModelling.jl/stable/#Process-API)
