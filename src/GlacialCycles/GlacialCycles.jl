module GlacialCycles
using ConceptualClimateModels

# TODO
for (root, dirs, files) in walkdir(joinpath(@__DIR__))
    file == "GlacialCycles.jl" && continue
    for file in files
        include(joinpath(root, file))
    end
end

end