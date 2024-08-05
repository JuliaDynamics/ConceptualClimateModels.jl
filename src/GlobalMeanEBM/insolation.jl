export AstronomicalForcingDeSaedeleer, astronomical_forcing_desaedeleer

"""
    AstronomicalForcingDeSaedeleer(; S = S, extensive = false)

Create the equation `S ~ astronomical_forcing_desaedeleer(t, extensive)`
which is Eq. (1) of [DeSaedeleer2013](@cite):

```@math
S = \\sum_i s_i \\sin(\\omega_i t) + c_i \\cos(\\omega_i t)
```

where the values of ``\\omega_i, s_i, c_i`` come from [Berger1978](@cite) who performed
a spectral expansion of the insolation. The validity range of this approximation is [-1, 0] Myr.

In the summation ``i`` goes up to 35 if `extensive`, otherwise up to 8.
The components are sorted according to magnitude of the spectral line,
so the default version has only the 8 most important spectra lines.

Note that in contrast to Eq. (1) of [DeSaedeleer2013](@cite) we do not normalize
``f`` and its value is in W/m² (the mean value is still deducted).
Additionally, the values of ``\\omega_i`` have been adjusted to expect
time in units of seconds.
"""
function AstronomicalForcingDeSaedeleer(; S = S, kw...)
    return S ~ astronomical_forcing_desaedeleer(t, kw...)
end

function astronomical_forcing_desaedeleer(t, extensive = false)
    F = 0.0
    limit = extensive ? 35 : 8
    @inbounds for i in 1:limit
        ω, s, c = @view(ORBITAL_FORCING_COMPONENTS[i, :])
        ωt = ω*t
        F += s*sin(ωt) + c*cos(ωt)
    end
    return F # do not scale, return it in W/m^2
end

# we need to register here otherwise the symbolic expression is too large,
# making everything too inefficient
@register_symbolic astronomical_forcing_desaedeleer(t)
@register_symbolic astronomical_forcing_desaedeleer(t, e::Bool) false


"""
    ORBITAL_FORCING_COMPONENTS

The full Table 1 of [DeSaedeleer2013](@ref) as a 35 by 3 matrix
with the column order the same as Table 1: ωᵢ sᵢ cᵢ and the difference
that column 1 is normalized to be in units of rad/sec instead of rad/kyr.
In contrast to the sorting of the original Table 1 here the elements are sorted by
amplitude of the spectral lines, i.e., sorted by `row[2]^2 + row[3]^2` in reverse
order, having highest amplitude (and hence most important) spectral lines first.
See the source code for the original sorting with obliquity terms first and
then precession terms.
"""
const ORBITAL_FORCING_COMPONENTS = begin
    ORIGINAL = [
    # Obliquity terms
    0.153249478547167 -11.2287376815124 3.51682075211241
    0.158148666238883 -3.82499371467540 -0.761851750263805
    0.117190147169570 2.28814805956066 1.80233702684623
    0.155061775112933 -1.29770081956440 -0.635152963728496
    0.217333905941751 0.380973541305497 -1.46301711999210
    0.150162587421217 1.54904176353302 -0.0883941912769817
    0.211709630908568 -0.810768209286259 -0.577980646565494
    0.156336369673117 -0.918358442095885 0.196083726889428
    0.148350290855451 0.256895610735773 -0.524697312305024
    0.206924898030688 -0.335783913402678 -0.0194792150128644
    0.212525165090383 0.267659228540196 0.128915417116900
    0.229992875969202 0.0696189733188958 0.0746231714061285
    0.306498957094334 0.0247349748169616 0.0140464395340974
    0.311398144786051 0.0138353727621181 0.0304736668840422
    0.004899187691716 -0.160479848721994 0.0594077968934257
    # Precession terms
    0.264933601588513 -15.5490493322904 -9.70406287110532
    0.280151350350945 15.4319556361701 4.75247271131525
    0.331110950251899 9.0992249352734 -10.6115244887390
    0.328024059125949 -7.87065384013669 6.61544246063503
    0.326211762560183 0.813786144754451 -4.52641408099246
    0.269742342439881 0.0690448504314857 -3.31639260969558
    0.332923246817665 1.44050770785967 1.06339286050120
    0.371638925683567 0.925324276580528 -1.02066758672154
    0.275366617473065 0.997628846513796 -0.362906496840039
    0.323124871434233 -0.378637986107629 0.527217891742183
    0.259396912994958 0.339477750517033 -0.560509461538342
    0.324937167999999 -0.576082669762308 1.18669572739338
    0.334197841377850 0.346906064369828 -0.648189701487285
    0.274551083291250 -0.441772417569753 0.289576210423804
    0.418183080135680 -0.0184884064645011 0.109632390175297
    0.111684123041346 -0.428006728186239 0.357006342316690
    0.433400828898112 -0.0049199219454561 -0.106148873639336
    0.126901871803777 0.257509918217341 -0.377639794223366
    0.336010137943616 -0.421809264016129 0.324327509437558
    0.177861471704732 -0.161827722328271 -0.362683869407858
    ]
    amplitudes = map(row -> sqrt(row[2]^2 + row[3]^2), eachrow(ORIGINAL))
    sorted = ORIGINAL[sortperm(amplitudes; rev = true), :]
    # convert frequencies to units of 1/sec
    sorted[:, 1] .*= 1000sec_in_year
    return sorted
end

# formula from North 1975 theory of energy balance models
function insolation_at_latitude(lat) # lat in degrees
    x = sind(lat)
    p2(x) = (3x^2 - 1)/2
    s = 1 -0.482p2(x)
    return s*solar_constant
end
