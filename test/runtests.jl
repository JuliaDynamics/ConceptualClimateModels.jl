using ConceptualClimateModels
using ConceptualClimateModels.CCMV
using Test

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

ds = processes_to_coupledodes(p1)
mtk = referrenced_sciml_model(ds)
eqs = all_equations(mtk)

# This also tests a bunch of default processes
for var in (T, ε, ε1, ε2, ε3, ε4, q)
    @test has_symbolic_var(eqs, var)
end
@test has_symbolic_var(eqs, ε)
@test has_symbolic_var(eqs, :ε1_T_tanh_ref)
@test has_symbolic_var(eqs, mtk.ε_0)
@test !has_symbolic_var(equations(mtk), α)

@variables begin
    OLR1(t) = 10.0
    OLR2(t) = 10.0
    OLR3(t) = 10.0
    a1(t) = 0.1
    a2(t) = 0.1
    a3(t) = 0.1
    a4(t) = 0.1
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
]

ds = processes_to_coupledodes(p2)
mtk = referrenced_sciml_model(ds)
for var in (T, C, OLR1, OLR2, :α_bg, a1, a2, a3, a4)
    if has_symbolic_var(mtk, var)
        @test true
    else
        error("var $var doesn't exist")
    end
end
