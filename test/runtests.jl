using ConceptualClimateModels
using ConceptualClimateModels.GlobalMeanEBM
using Test

@testset "all processes" begin
function test_symbolic_var(mtk, var)
    if has_symbolic_var(mtk, var)
        @test true
    else
        error("var $var doesn't exist")
    end
end

@variables begin
    ε1(t) = 0.1
    ε2(t) = 0.1
    ε3(t) = 0.1
    ε4(t) = 0.1
end

p1 = [
    BasicRadiationBalance(),
    EmissivityStefanBoltzmanOLR(),
    IceAlbedoFeedback(; min = 0.3, max = 0.7),
    α ~ α_ice,
    ParameterProcess(ε), # emissivity is a parameter
    f ~ 0, # no external forcing
    # absorbed solar radiation has a default process
    EmissivityFeedbackTanh(ε = ε1),
    SoedergrenClearSkyEmissivity(ε = ε2),
    EmissivitySellers1969(ε = ε3),
    EmissivityLogConcentration(ε = ε4),
]

ds = processes_to_coupledodes(p1, GlobalMeanEBM)
mtk = referrenced_sciml_model(ds)
eqs = all_equations(mtk)

# This also tests a bunch of default processes
for var in (T, ε, ε1, ε2, ε3, ε4, q, :ε1_T_sigmoid_ref, mtk.ε_0)
    test_symbolic_var(eqs, var)
end


@variables begin
    OLR1(t) = 10.0
    OLR2(t) = 10.0
    OLR3(t) = 10.0
    a1(t) = 0.1
    a2(t) = 0.1
    a3(t) = 0.1
    a4(t) = 0.1
    ΔS1(t) = 0.0
    ΔT1(t) = 0.0
end

p2 = [
    LinearOLR(),
    LinearClearSkyOLR(; OLR = OLR1),
    BudykoOLR(; OLR = OLR2),
    CloudAlbedoExponential(),
    CloudAlbedoLinear(; α_cloud = a1),
    DirectAlbedoAddition(α = a4),
    CoAlbedoProduct(α = a2),
    SeparatedClearAllSkyAlbedo(α = a3),
    C ~ 0.6,
    ΔTLinearRelaxation(),
    ΔTStommelModel(ΔT = ΔT1, ΔS = ΔS1),
]

ds = processes_to_coupledodes(p2, GlobalMeanEBM)
mtk = referrenced_sciml_model(ds)
for var in (T, C, OLR1, OLR2, :α_bg, a1, a2, a3, a4, ΔS1, ΔT)
    if has_symbolic_var(mtk, var)
        @test true
    else
        error("var $var doesn't exist")
    end
end



end


@testset "dynamical systems integration" begin
    using DynamicalSystems

    budyko_processes = [
        BasicRadiationBalance(),
        EmissivityStefanBoltzmanOLR(),
        IceAlbedoFeedback(; min = 0.3, max = 0.7),
        α ~ α_ice,
        ParameterProcess(ε), # emissivity is a parameter
        f ~ 0, # no external forcing
        # absorbed solar radiation has a default process
    ]

    budyko = processes_to_coupledodes(budyko_processes, GlobalMeanEBM)

    grid = plausible_grid(budyko)
    mapper = AttractorsViaRecurrences(budyko, grid)
    rfam = RecurrencesFindAndMatch(mapper)
    sampler = plausible_ic_sampler(budyko)

    set_parameter!(budyko, :ε_0, 0.5)

    fracs = basins_fractions(mapper, sampler)
    attractors = extract_attractors(mapper)

    @test length(attractors) == 2

end