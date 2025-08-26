# [CloudToppedMixedLayerModel](@id ctmlm_examples)

```@example MAIN
using DynamicalSystems, ConceptualClimateModels
import ConceptualClimateModels.CloudToppedMixedLayerModel as CTMLM
```

## Fixed states of Stevens 2006

Here we want to recreate the same bulk boundary layer model defined in [Stevens2006](@cite),
equations 31-33.

```@example MAIN
const cₚ = 1004.0 # heat capacity at constant pressure (J/K/kg)

eqs = [
    CTMLM.mlm_dynamic(),
    CTMLM.entrainment_velocity(:Stevens2006; use_augmentation = false),
    ## Figure 1 and Sec. 4.2 provide these values:
    CTMLM.s₊ ~ 301200.0/cₚ,
    CTMLM.q₊ ~ 1.56,
    CTMLM.s₀ ~ CTMLM.s₊ - 12.5,
    CTMLM.ρ₀ ~ 1,
    CTMLM.q₀ ~ CTMLM.q_saturation(288.96 + 1.25),
    CTMLM.ΔF_s ~ 40.0
]
```

We now give these equations to the main function that creates a `DynamicalSystem` from the processes (and we also provide `CTMLM` to obtain default processes from)

```@example MAIN
ds = processes_to_coupledodes(eqs, CTMLM)
```

We also make sure to use the same parameter values as in the paper:
```@example MAIN
set_parameter!(ds, :D, 4e-6)
set_parameter!(ds, :c_d, 0.0011)
set_parameter!(ds, :U, 0.008/0.0011) # U*c_d = 0.008
set_parameter!(ds, :e_e, 1.0) # effective emissivity
```

And that's all. We can run the system to its steady state:

```@example MAIN
step!(ds, 100.0)
```

leading to a steady state with height ``z_b`` (``h_\infty`` in [Stevens2006](@cite)) of about 800 meters as in the paper. Extracting the variable ``\sigma`` of Eq. 38 is also very easy:

```@example MAIN
Δs = CTMLM.s₊ - CTMLM.s₀ # symbolic variable not existing in the graph of the `ds`
σ = observe_state(ds, CTMLM.V * Δs / CTMLM.ΔF_s * cₚ)
```
and we get the same value (note multiplication by `cₚ`, because `s` is in units of K).

If you want to see a list of all equations that compose the dynamical system then do

```@example MAIN
all_equations(ds)
```

In the compiled documentation these render via LaTeX but running this in the REPL won't be as pretty :)

## Adding clouds

We can add clouds to the mixed layer model using the same decoupling-based approach as in [Datseris2025](@cite) (while keeping SST a fixed boundary condition)
just by including a couple more equations to the ones already defined,
so that ``C, \mathcal{D}`` have a process assigned to them.
We will also augment `ΔF` to be partially proportional to `C` using a very simple ad-hoc approach.

```@example MAIN
eqs = [
    CTMLM.mlm_dynamic(),
    CTMLM.entrainment_velocity(:Stevens2006; use_augmentation = false),
    ## Cloud stuff
    CTMLM.cf_dynamic(),
    CTMLM.decoupling_variable(),
    CTMLM.cloud_base_height(:Bolton1980),
    CTMLM.CTRC ~ 10 + 40*CTMLM.C, # cloud top radiative cooling
    CTMLM.ΔF_s ~ CTMLM.CTRC, # same equation
    ## rest the same
    CTMLM.s₊ ~ 301200.0/cₚ,
    CTMLM.q₊ ~ 1.56,
    CTMLM.s₀ ~ CTMLM.s₊ - 12.5,
    CTMLM.ρ₀ ~ 1,
    CTMLM.q₀ ~ CTMLM.q_saturation(288.96 + 1.25),
]

ds = processes_to_coupledodes(eqs, CTMLM)
```

And we run the model to a steady state:
```@example MAIN
set_parameter!(ds, :D, 4e-6)
set_parameter!(ds, :c_d, 0.0009)
set_parameter!(ds, :U, 6.8) # U*c_d = 0.008
set_parameter!(ds, :e_e, 1.0) # effective emissivity

step!(ds, 100.0)
observe_state.(ds, (:z_b, :q_b, :s_b, :C)) # get state variables by name
```

## Climate change scenario: increasing SST

If you want to study multistability for alternate cloud states (Cumulus vs Stratocumulus), or perform continuations (like the climate change scenarios in [Datseris2025](@cite)), visit the documentation of Attractors.jl.

As a brief example we will perform a simple climate change scenario where SST increases with a constant rate. First, we modify the equations so that `s₀, q₀` are derived from a prescribed SST:

```@example MAIN
@parameters SST = 290.0
eqs = [
    CTMLM.mlm_dynamic(),
    CTMLM.entrainment_velocity(:Stevens2006; use_augmentation = false),
    ## Cloud stuff
    CTMLM.cf_dynamic(),
    CTMLM.decoupling_variable(),
    CTMLM.cloud_base_height(:Bolton1980),
    CTMLM.CTRC ~ 10 + 40*CTMLM.C, # cloud top radiative cooling
    CTMLM.ΔF_s ~ CTMLM.CTRC, # same equation
    ## Usage of SST
    CTMLM.s₀ ~ SST,
    CTMLM.q₀ ~ CTMLM.q_saturation(SST),
    ## rest the same
    CTMLM.s₊ ~ 301200.0/cₚ,
    CTMLM.q₊ ~ 1.56,
    CTMLM.ρ₀ ~ 1,
]

ds = processes_to_coupledodes(eqs, CTMLM)
set_parameter!(ds, :D, 4e-6)
set_parameter!(ds, :c_d, 0.0009)
set_parameter!(ds, :U, 6.8) # U*c_d = 0.008
set_parameter!(ds, :e_e, 1.0) # effective emissivity
ds
```

To perform global continuation we need to create an attractor mapper, which we create here by tessellating the state space. Have a look at the Attractors.jl tutorial if the next few lines of code are puzzling to you

```@example MAIN

grid = (
    (1:1.0:25), # specific humidity
    (270:5.0:330), # static energy
    (0:100.0:3000), # height
    (0:0.1:1), # cloud fraction
)

sampler, = statespace_sampler(grid)

mapper = AttractorsViaRecurrences(ds, grid)
ascm = AttractorSeedContinueMatch(mapper)

fs = basins_fractions(mapper, sampler; show_progress = false)
attractors = extract_attractors(mapper)
```

with the `sampler, mapper, ascm` data structures in order, we can easily now run a global continuation with changing SST:

```@example MAIN
prange = 290:1:310
pidx = :SST
fractions_cont, attractors_cont = global_continuation(ascm, prange, pidx, sampler; show_progress = false)
attractors_cont
```

if you are not sure what the output means, no worries, just have a look at the Attractors.jl documentation. Here we visualize the cloud fraction:

```@example MAIN
using CairoMakie

# cloud fraction and height values of last point on the attractor
varvalues = [A -> A[end][4], A -> A[end][3]]

fig = plot_basins_attractors_curves(
	fractions_cont, attractors_cont, varvalues, prange,
)
content(fig[2,1]).ylabel = "C"
content(fig[3,1]).ylabel = "z_b"
fig
```

we see for our ad hoc parameterizations, the dynamical system has no stable states for SST > 297. Before that cloud fraction decreases very fast from a Stratocumulus state to a Cumulus one.

The model diverges due to the definition of the entrainment velocity (which is proportional to `ΔF_s/Δ₊s`) in combination with a fixed `s₊`: the model reaches a state where `Δ₊s` is so small that height increases unnaturally well beyond the state space tessellation we defined.