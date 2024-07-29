module AdaptivePredicates

include("init.jl")
include("caches.jl")
include("macros.jl")
include("arithmetic.jl")
include("predicates.jl")

export orient2, orient3, incircle, insphere
export orient2p, orient3p, incirclep, inspherep
export orient2fast, orient3fast, incirclefast, inspherefast

end # module
