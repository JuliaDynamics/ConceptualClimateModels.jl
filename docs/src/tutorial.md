# [Tutorial](@id tutorial)

## Terminology

ConceptualClimateModels.jl follows a process-based modelling approach
to make differential equation systems from _processes_.
A _process_ is simply a particular _equation_ defining the dynamics of a climate
variable, while also explicitly defining _which_ variable the equation defines.
A vector of processes is composed by the user, and given to the main function [`processes_to_coupledodes`](@ref) which bundles them into a system of equations
that creates a dynamical system. The dynamical system can then be further analyzed
in terms of stability properties, multistability, tipping, periodic (or not)
behavior, and many more aspects, via the DynamicalSystems.jl library (see the examples).

Note the distinction: a _process_ is _not_ the climate variable
(such as "clouds" or  "insolation"); rather it is the _exact equation_ that
defines the behavior of the climate variable, and could itself utilize
many other already existing climate variables, or introduce new ones.
In terminology of climate science a _process_ is a generalization
of the term "parameterization".
Many different processes may _describe_ the behavior of a particular variable and
typically one wants to analyze what the model does when using one versus the other
process.


!!! note "Familiarity with DynamicalSystems.jl and ModelingToolkit.jl"
    ConceptualClimateModels.jl builds on [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/) for building the equations representing the climate model, and it builds on [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/dev/) to further analyze the models. Familiarity with either package is good to have, and will allow you to faster and better understand the concepts discussed here. Nevertheless familiarity is actually optional as the steps required to use ConceptualClimateModels.jl are all explained in this tutorial.


## Introductory example

Let's say that we want to create the most basic energy balance model,
```math
c_T \frac{dT}{dt} = ASR - OLR + f
```
where ``ASR`` is the absorbed solar radiation given by
```math
ASR = S (1-\alpha)
```
with ``\alpha`` the planetary albedo, ``OLR`` is the outgoing longwave radiation given by the linearized expression
```math
OLR = A + BT
```
and ``f`` some radiative forcing at the top of the atmosphere, that is based on CO2 concentrations and given by
```math
f = 3.7\log_2\left(\frac{CO_2}{400}\right)
```
with CO2 concentrations in ppm.

To create this model with ConceptualClimateModels.jl while providing the least information possible we can do:

```@example MAIN
using ConceptualClimateModels
using ConceptualClimateModels.GlobalMeanEBM # bring names into scope

processes = [
    BasicRadiationBalance(),
    LinearOLR(),
    ParameterProcess(α),
    CO2Forcing(), # note that for default CO2 value this is zero forcing
]

ds = processes_to_coupledodes(processes, GlobalMeanEBM)
println(dynamical_system_summary(ds))
```

The output is a dynamical system from [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/stable/) that is generated via symbolic expressions based on [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/), utilizing the process-based approach of [ProcessBasedModelling.jl](https://juliadynamics.github.io/ProcessBasedModelling.jl/stable/).

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

The symbolic variables and parameters can also be used to query or alter the
dynamical system. For example, we can obtain, or alter, any parameter by providing the
symbolic parameter index.
There are multiple ways to obtain the symbolic index provided we know its name.
First, we can re-create the symbolic parameter:

```@example MAIN
index = first(@parameters CO2_0)
```

Second, we can use the retrieved MTK model and access its `CO2_0` parameter:


```@example MAIN
index = mtk.CO2_0
```

Third, we can use a `Symbol` corresponding to the variable name.
This is typically the simplest way.
```@example MAIN
index = :CO2_0
```

we can query the value of this named parameter in the system,

```@example MAIN
current_parameter(ds, index)
```

or alter it:

```@example MAIN
# access symbolic parameter CO2_0 from the tracked symbolic list of the model
set_parameter!(ds, index, 800)
```

Similarly, we can obtain or alter values corresponding to the dynamic variables,
or observed functions of the state of the system, using their symbolic indices.
For example we can obtain the value corresponding to symbolic variable ``T`` by:

```@example MAIN
observe_state(ds, T) # binding `T` already exists in scope
```

or obtain the ``OLR`` (outgoing longwave radiation)

```@example MAIN
observe_state(ds, :OLR)
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
using ProcessBasedModelling: lhs, rhs
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

# Default processes

Hold on a minute though, because in the original processes we provided,
```julia
processes = [
    BasicRadiationBalance(),
    LinearOLR(),
    ParameterProcess(α),
    CO2Forcing(),
]
```
there was no process that defined for example the absorbed solar radiation ``ASR``!

Well, ConceptualClimateModels.jl allows the concept of "default processes".
The package exports some **submodules**, and each submodule targets a particular application of conceptual climate models. In this Tutorial we are using the [`GlobalMeanEBM`](@ref) submodule, which provides functionality to model global mean energy balance models.

Each submodule defines and exports its own **list of predefined symbolic variables**.
When we wrote
```julia
using ConceptualClimateModels.GlobalMeanEBM
```
we brought into scope all the variables that this (sub)module defines and exports, such as `T, α, OLR`.

Each (sub)module and provides a list of **predefined processes for its predefined symbolic variables**.
These predefined default processes are loaded automatically when we provide the (sub)module as a second argument to `processes_to_coupledodes`, which we did above.
In this (sub)module, the default process for the $ASR$ is $ASR = S(1-\alpha)$ with $S$ the solar constant.

The function [`processes_to_coupledodes`](@ref) goes through all processes the user provided and identifies variables that themselves do not have a process.
It then checks the list of default processes and attempts to assign one to these variables.
If there are no default processes for some variables, it makes the variables themselves parameters with the same name but with a subscript 0 if the variables
have a default value.

For example, let's assume that we completely remove default processes and we don't specify
a process for the albedo ``α``:

```@example MAIN
processes = [
    BasicRadiationBalance(),
    LinearOLR(),
    CO2Forcing(), # note that for default CO2 values this is zero forcing
    ASR ~ S*(1-α), # add the processes for ASR, but not for S or α
]

# note the empty list as 2nd argument, which is otherwise
# the default processes. Notice also that we make an MTK model
# (which is the step _before_ converting to a dynamical system)
mtk = processes_to_mtkmodel(processes, [])
# we access the equations directly from the model
equations(mtk)
```

You will notice the equation ``α = α_0``
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
└ `ParameterProcess(α)`.
```

[`ParameterProcess`](@ref) is the most trivial process: it simply means that the corresponding variable does not have any physical process and rather it is a system
parameter.

This automation does not occur if there is no default value.
For example, the $ASR$ variable does not have a default value. If we have not assigned a process for $ASR$, the system construction would error instead:

```julia
processes = [
    BasicRadiationBalance(),
    LinearOLR(),
    CO2Forcing(),
]

mtk = processes_to_mtkmodel(processes, [])
```

```
ERROR: ArgumentError: Variable ASR(t) was introduced in process of
variable T(t). However, a process for ASR(t) was not provided,
there is no default process for ASR(t), and ASR(t) doesn't have a default value.
Please provide a process for variable ASR(t).
```

These warnings and errors are always "perfectly" informative.
They tell us exactly which variable does not have a process, and exactly which other process introduced the process-less variable first.
This makes the modelling experience stress-free, especially when large and complex
models are being created.

## Adding your own processes

Each of the submodules of ConceptualClimateModels.jl provides an increasing list of
predefined processes that you can use out of
the box to compose climate models. The predefined processes all come from existing
literature and cite their source via BiBTeX.

It is also very easy to make new processes on your own.
The simplest way to make a process is to just provide an equation for it
with the l.h.s. of the equation being the variable the process defines.

For example:

```@example MAIN
@variables x(t) = 0.5 # all variables must be functions of (t)
x_process = x ~ 0.5*T^2 # x is just a function of temperature
```

A more re-usable approach however is to create a function that generates a process
or create a new process type as we describe in [making new processes](@ref new_processes).

## A note on symbolic variable instances

Recall that when we wrote
```julia
using ConceptualClimateModels.GlobalMeanEBM
```
at the start of this tutorial, we brought into scope variables that this (sub)module defines and exports, such as `T, α, OLR`. They are listed on the (sub)module's documentation page, and are used in that module's default processes.

In the (sub)module predefined processes you will notice call signatures like
```julia
BasicRadiationBalance(; T, f, kwargs...)
```
There are keywords that do not have an assignment like `T, f` above.
They use the [(sub)module's predefiend variables.

Crucially, these default variables are _symbolic variables_. They are defined as
```julia
@variables begin
    T(t) = 0.5 # ...
    # ...
end
```
which means that expressions that involve them result in symbolic expressions,
for example
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

You can use your own variables in any predefined processes.
You can define them by doing
```@example MAIN
@variables begin
    (T1_tropics(t) = 290.0), [bounds = (200.0, 350.0), description = "temperature in tropical box 1, in Kelvin"]
end
```
and then assign them to the corresponding keyword argument
```@example MAIN
process = BasicRadiationBalance(T = T1_tropics)
```

Defining variables with the extra `bounds, description` annotations is
useful for integrating with the rest of the functionality of the library,
and therefore it is strongly recommended.

!!! warn "Custom variables need to be assigned everywhere!"
    Above we assigned `T1_tropics` as the temperature variable.
    This means we also need to assign the same variable as the one setting the
    `OLR` variable by also providing the processes
    `LinearOLR(T = T1_tropics)` (for example).


## Default values, limits, etc.

All predefined variables that could be dynamic variables (i.e., could have a time derivative applied to them) have a default value, a description, and plausible physical bounds.

To obtain the default value, use `default_value(x)`. For the description,
use `getdescription(x)`. For the bounds, use:

```@docs
physically_plausible_limits(::Any)
```

e.g.,

```@example MAIN
physically_plausible_limits(T)
```

## Automatically named parameters

The majority of predefined processes create symbolic
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
    α_ice = my_ice_α,
    Tfreeze = myfreeze, # my custom parameter
    max = LiteralParameter(0.9) # don't make a parameter
)

mtk = processes_to_mtkmodel([ice_process])
equations(mtk)
```

```@example MAIN
parameters(mtk)
```

## Integration with DynamicalSystems.jl

ConceptualClimateModels.jl integrates with [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/dev/) by providing initial condition
sampling to use when e.g., finding attractors and their basin fractions with
`DynamicalSystems.basins_fractions`, and with the function [`dynamical_system_summary`](@ref).
Moreover, since all dynamical systems generated by ConceptualClimateModels.jl have
symbolic bindings, one can use the symbolic variables in e.g., [interactive GUI exploration](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/dev/visualizations/#Interactive-or-animated-trajectory-evolution)
or to access or set the parameters of the system.

```@docs
physically_plausible_limits(::DynamicalSystem)
physically_plausible_ic_sampler
plausible_grid
dynamical_system_summary
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
   or read the documentation string of [`Process`](@ref) below.

```@docs
Process
```
