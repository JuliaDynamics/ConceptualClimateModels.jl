export EmissivityStefanBoltzmanOLR
export LinearOLR, LinearClearSkyOLR

"""
    EmissivityStefanBoltzmanOLR(; ε = ε, T = T)

Create the equation `OLR ~ ε*σ*T^4` where `σ` is the Stefan Boltzmann constant
and `ε` the effective emissivity, also known as the "grayness" of the system,
or the deviation it has from being a perfect black body [Ghil1981](@cite).
`ε` then needs to be parameterized itself to include greenhouse or other climate effects.
"""
function EmissivityStefanBoltzmanOLR(; ε = ε, T = T)
    return OLR ~ ε*(σ_Stefan_Boltzman)*T^4
end

"""
    LinearOLR(; T = T, A = -277.0, B = 1.8)

Create the equation `OLR ~ A + B*T`.
This is a linearized outgoing longwave radiation (OLR), and is the
same equation as (7) of [North1981](@cite): ``OLR = A + BT``
with ``T`` temperature in Kelvin and ``A, B`` constants.
However, default ``A, B`` are fitted from current CERES all sky OLR
and using ERA5 data for the 2-meter temperature. We assume `T` in Kelvin.
This linear approximation is quite accurate for temporally averaged data
``T \\in (220, 280)`` however drops drastically in accuracy after that
due to the nonlinear effects of clouds (as evident by observational data).

We note a big difference between current CERES data and the values reported
in [North1981](@cite): here `A=214.67` (assuming ``T`` in Celcius) and `B=1.8`
versus the values `A=203.3` and `B=2.09` in [North1981](@cite).

If instead of all sky, if we fit the clear sky CERES data, we get
`A = -326.0, B = 2.09`.
Interestingly, coefficient `B` here is the same as that reported by [North1981](@cite),
but `A=244.88` (assuming `T` in Celcius) is not.
"""
function LinearOLR(; T = T, A = -277.0, B = 1.8)
    return OLR ~ A + B*T
end


"""
    LinearClearSkyOLR(; T = T)

Equivalent with `LinearOLR(; T, A = A = -326.0, B = 2.09)` and provided
as a convenience for the clear sky fit to CERES data.
"""
function LinearClearSkyOLR(; T = T)
    return LinearOLR(; T, A = -326.0, B = -2.09)
end
