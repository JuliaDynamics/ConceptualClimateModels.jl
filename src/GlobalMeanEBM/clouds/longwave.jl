export BudykoOLR

"""
    BudykoOLR(; OLR=OLR, T=T, C=C,
        BudykoOLR_A = -461.8068, BudykoOLR_B = 2.58978,
        BudykoOLR_Ac = -377.22741, BudykoOLR_Bc = 1.536171
    )

Create the equation `OLR ~ A + B*T - C*(Ac + Bc*T)` for the dependence of OLR on
both temperature and cloud fraction (in 0-1). This is the same
as Eq. (1) of [Budyko1969](@cite). However, here `T` is expected in Kelvin,
and the coefficients have been extracted by fitting into CERES data in the same
way as in [`LinearOLR`](@ref).
"""
function BudykoOLR(; OLR=OLR, T=T, C=C,
        BudykoOLR_A = -461.8068, BudykoOLR_B = 2.58978,
        BudykoOLR_Ac = -377.22741, BudykoOLR_Bc = 1.536171
    )
    @convert_to_parameters BudykoOLR_A BudykoOLR_B BudykoOLR_Ac BudykoOLR_Bc
    return OLR ~ BudykoOLR_A + BudykoOLR_B*T - C*(BudykoOLR_Ac + BudykoOLR_Bc*T)
end