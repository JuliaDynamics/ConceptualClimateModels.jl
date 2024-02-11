abstract type AlbedoProcess <: Process end

export DirectAlbedoAddition, CoAlbedoProduct, SeparatedClearAllSkyAlbedo

ProcessBasedModelling.lhs_variable(::AlbedoProcess) = α

"""
    DirectAlbedoAddition(; α_bg = 0.1, other_albedo_variables = (α_ice, α_clouds))

Create the equation `α ~ α_bg + other_albedo_variables...`, meaning that
planetary albedo `α` is a direct sum of all specified albedos.
"""
Base.@kwdef struct DirectAlbedoAddition <: AlbedoProcess
    α_bg::Real = 0.1
    other_albedo_variables = (α_ice, α_cloud)
end
function ProcessBasedModelling.rhs(p::DirectAlbedoAddition)
    @parameters α_bg = p.α_bg
    # Add more other_albedo_variables here if any exist
    return α_bg + (p.other_albedo_variables...)
end

"""
    CoAlbedoProduct(albedo_variables = (α_ice, α_cloud))

Create the equation `1 - α ~ prod(a -> (1 - a), albedo_variables)`
meaning that the planetary co-albedo is the product of the co-albedos of all
albedo variables. This would be the planetary albedo
if all components were uniform layers, while the bottom-most layer (surface)
had perfect absorption and all other layers had 0 absorption and finite reflection.
"""
Base.@kwdef struct CoAlbedoProduct <: AlbedoProcess
    albedo_variables = (α_ice, α_cloud)
end
function ProcessBasedModelling.rhs(p::CoAlbedoProduct)
    a = p.albedo_variables
    return 1 - prod(α -> (1.0 - α), a)
end

"""
    SeparatedClearAllSkyAlbedo(; α_cloud = α_cloud, C = C, α_clr = 0.15)

Create the equation `α ~ α_cloud*C + α_clr*(1 - C)`.

[Bender2017](@cite) argue that one can assume a separation between clear-sky
and cloud albedo, so that `α = α_cloud*C + α_clr*(1 - C)` with `C` the cloud fraction
and `α_clr` the clear sky albedo. They further cite [Cess1976](@cite) to facilitate the claim
Additionally, Eq. (20) of [Barker1999](@cite) provides an identical expression.
"""
function SeparatedClearAllSkyAlbedo(; α_cloud = α_cloud, C = C, α_clr = 0.15)
    return α ~ α_cloud*C + α_clr*(1 - C)
end