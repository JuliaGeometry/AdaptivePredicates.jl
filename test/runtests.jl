using Test
import ExactPredicates
using Supposition
using AdaptivePredicates

cd("original") do
    if !ispath("libpredicates.so")
        run(`gcc -shared -fPIC -O3 predicates.c -o libpredicates.so`)
    end
end
include("original/CPredicates.jl")

function orient(a::Complex, b::Complex, c::Complex)
    orient((a.re, a.im), (b.re, b.im), (c.re, c.im))
end

function orient(a, b, c)
    dir = orient2(a, b, c)
    if dir < 0.0
        return -1
    elseif dir > 0.0
        return 1
    else
        return 0
    end
end

@testset "Spot checks found through Supposition" begin
    (a,b,c) = (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)
    @test orient(a,b,c) == ExactPredicates.orient(a,b,c)
    @test orient(a,b,c) == CPredicates.orient(a, b, c)

    (a, b, c) = (-0.08908073736089261 - 0.1026469210893071im, -1.956194437297279 - 2.254105006193488im, -0.3200306149494339 - 0.36876835836900623im)
    @test orient(a,b,c) == CPredicates.orient(a, b, c)
    @test orient(a,b,c) == ExactPredicates.orient(a,b,c)

    (a, b, c) = (a = 0.0 + 0.0im, b = 0.0 + 5.0e-324im, c = 5.0e-324 + 2.2536010670724063e30im)
    @test orient(a,b,c) == CPredicates.orient(a, b, c)
    @test_broken orient(a,b,c) == ExactPredicates.orient(a,b,c)

    (a,b,c) = (a = 0.0 + 0.0im, b = 0.0 + 5.0e-324im, c = 5.060829986287029e79 + 3.552170572284382e228im)
    @test orient(a,b,c) == CPredicates.orient(a, b, c)
    @test_broken orient(a,b,c) == ExactPredicates.orient(a,b,c)

    (a, b, c) = (a = 0.0 + 0.0im, b = 0.0 + 5.0e-324im, c = 5.0e-324 + 0.0im)
    @test orient(a,b,c) == CPredicates.orient(a, b, c)
    @test_broken orient(a,b,c) == ExactPredicates.orient(a,b,c)

    (a,b,c) = (a = 0.0 + 0.0im, b = 0.0 + 5.495397811658139e-246im, c = -4.495267265324774e-79 - 1.382478711138964e-74im)
    @test orient(a,b,c) == CPredicates.orient(a, b, c)
    @test orient(a,b,c) == ExactPredicates.orient(a,b,c)
end

const floatgen = Data.Floats{Float64}()
const complexgen = @composed function _complex(a=floatgen, b=floatgen)
    a+b*im
end

# @check function orient_c(a=complexgen, b=complexgen, c=complexgen)
#     dir = CPredicates.orient(a,b,c)
#     event!("Adaptive", dir)
#     edir = CPredicates.orient_exact(a,b,c)
#     event!("Exact", edir)
#     dir == edir
# end

@check function output(a=complexgen, b=complexgen, c=complexgen)
    orient(a,b,c) ∈ (-1, 0, 1)
end

# Note: This should be true for all, if not this is a bug in the port
@check function orient_against_c(a=complexgen, b=complexgen, c=complexgen)
    dir = orient(a,b,c)
    event!("Port", dir)
    cdir = CPredicates.orient(a,b,c)
    event!("Original", cdir)
    dir == cdir
end

# The next two should agree
@check function orient_against_exact(a=complexgen, b=complexgen, c=complexgen)
    dir = orient(a,b,c)
    event!("Adaptive", dir)
    edir =  ExactPredicates.orient(a,b,c)
    event!("Exact", edir)
    dir == edir
end

@check function orient_c_against_exact(a=complexgen, b=complexgen, c=complexgen)
    dir = CPredicates.orient(a,b,c)
    event!("Original", dir)
    edir =  ExactPredicates.orient(a,b,c)
    event!("Exact", edir)
    dir == edir
end

