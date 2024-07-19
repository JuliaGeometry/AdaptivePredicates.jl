module AdaptivePredicates

include("init.jl")
include("caches.jl")
include("macros.jl")
include("arithmetic.jl")
include("predicates.jl")

export orient2, orient3, incircle, insphere
export orient2p, orient3p, incirclep, inspherep

@static if VERSION ≥ v"1.11.0-DEV.469"
    eval(Meta.parse("public orient2, orient3, incircle, insphere, orient2p, orient3p, incirclep, inspherep"))
end

end # module