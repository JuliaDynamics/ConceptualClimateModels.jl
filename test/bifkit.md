### Using BifurcationKit.jl

Because for this example having _unstable_
fixed point branches is crucial, we will perform the alternative continuation
provided by BifurcationKit.jl. BifurcationKit.jl does not provide most of the conveniences
that DynamicalSystems.jl does. For example, it does not integrate well enough with
DifferentialEquations.jl (to allow passing `ODEProblem` which is created by `DynamicalSystem`).
It also does not allow indexing parameters by their symbolic bindings.

Due to this, we need to change the dynamical system generation to disable [parameter splitting](https://github.com/SciML/ModelingToolkit.jl/issues/2482#issuecomment-1959330130),
so that the parameter container becomes a simple vector of numbers
that can be indexed by the integers:

```@example MAIN
stommel = processes_to_coupledodes(processes; split = false, inplace = true)
```

We now need to find out which integer corresponds the the parameters of interest,
which we do by seeing the list:

```@example MAIN
mtk = referrenced_sciml_model(stommel)
pnames = parameters(mtk)
hcat(eachindex(pnames), pnames)
```

from where we obtain that index 1 corresponds to `r_η`.

From here, we follow the [introductory tutorial](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/gettingstarted/) of BifurcationKit.jl
and first make a `BifurcationProblem`:

```@example MAIN
import BifurcationKit

vf = dynamic_rule(stommel) # vector field, the dynamic rule
u = current_state(stommel) # initial guess for bifurcation
p = current_parameters(stommel) # parameter container
i = (BifurcationKit.@lens _[1]) # BifurcationKit.jl doesn't allow normal indices, only some weird "lenses"

prob = BifurcationKit.BifurcationProblem(vf, u, p, i)
```

Great, that took a while but now we can finally perform the continuation

```@example MAIN
opts = BifurcationKit.ContinuationPar(p_min = 0.0, p_max = 5.0)
br = BifurcationKit.continuation(prob, BifurcationKit.PALC(), opts)
```