export BudykoOLR

"""
    BudykoOLR(; T=T, C=C,
        A = -461.8068, B = 2.58978, Ac = -377.22741, Bc = 1.536171
    )

Create the equation `OLR ~ A + B*T - C*(Ac + Bc*T)` for the dependence of OLR on
both temperature and cloud fraction (in 0-1). This is the same
as Eq. (1) of [Budyko1969](@cite). However, here `T` is expected in Kelvin,
and the coefficients have been extracted by fitting into CERES data in the same
way as in [`LinearOLR`](@ref).
"""
function BudykoOLR(; T=T, C=C,
        A = -461.8068, B = 2.58978, Ac = -377.22741, Bc = 1.536171
    )
    @parameters begin
        Budyko_A = A
        Budyko_B = B
        Budyko_Ac = Ac
        Budyko_Bc = Bc
    end
    return OLR ~ Budyko_A + Budyko_B*T - C*(Budyko_Ac + Budyko_Bc*T)
end