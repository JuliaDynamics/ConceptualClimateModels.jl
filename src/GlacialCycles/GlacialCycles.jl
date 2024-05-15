module GlacialCycles
using ConceptualClimateModels

# TODO
for (root, dirs, files) in walkdir(joinpath(@__DIR__))
    for file in files
        file == "GlacialCycles.jl" && continue
        include(joinpath(root, file))
    end
end

end