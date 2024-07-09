using Test
import ExactPredicates
using Supposition
using AdaptivePredicates

cd("original") do
    if !isfile("libpredicates.so")
        if Sys.iswindows()
            # random(), being a POSIX function, is not available on Windows. We need to use rand().
            pred = read("predicates.c", String)
            new_pred = replace(pred, "random()" => "rand()")
            write("_predicates.c", new_pred)
            pred_path = "_predicates.c"
        else
            pred_path = "predicates.c"
        end
        run(`gcc -shared -fPIC -g2 $pred_path -o libpredicates.so`)
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
    (a, b, c) = (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)
    @test orient(a, b, c) == ExactPredicates.orient(a, b, c)
    @test orient(a, b, c) == CPredicates.orient(a, b, c)

    (a, b, c) = (-0.08908073736089261 - 0.1026469210893071im, -1.956194437297279 - 2.254105006193488im, -0.3200306149494339 - 0.36876835836900623im)
    @test orient(a, b, c) == CPredicates.orient(a, b, c)
    @test orient(a, b, c) == ExactPredicates.orient(a, b, c)

    (a, b, c) = (a=0.0 + 0.0im, b=0.0 + 5.0e-324im, c=5.0e-324 + 2.2536010670724063e30im)
    @test orient(a, b, c) == CPredicates.orient(a, b, c)
    @test_broken orient(a, b, c) == ExactPredicates.orient(a, b, c)

    (a, b, c) = (a=0.0 + 0.0im, b=0.0 + 5.0e-324im, c=5.060829986287029e79 + 3.552170572284382e228im)
    @test orient(a, b, c) == CPredicates.orient(a, b, c)
    @test_broken orient(a, b, c) == ExactPredicates.orient(a, b, c)

    (a, b, c) = (a=0.0 + 0.0im, b=0.0 + 5.0e-324im, c=5.0e-324 + 0.0im)
    @test orient(a, b, c) == CPredicates.orient(a, b, c)
    @test_broken orient(a, b, c) == ExactPredicates.orient(a, b, c)

    (a, b, c) = (a=0.0 + 0.0im, b=0.0 + 5.495397811658139e-246im, c=-4.495267265324774e-79 - 1.382478711138964e-74im)
    @test orient(a, b, c) == CPredicates.orient(a, b, c)
    @test orient(a, b, c) == ExactPredicates.orient(a, b, c)
end

const floatgen = Data.Floats{Float64}()
const floatgen2 = Data.Floats{Float64}(infs=false, nans=false) # IntervalArithmetic.jl screams out warnings with bad intervals during testing
const complexgen = @composed function _complex(a=floatgen, b=floatgen)
    a + b * im
end
const complexgen2 = @composed function _complex(a=floatgen2, b=floatgen2)
    a + b * im
end

@check function output(a=complexgen, b=complexgen, c=complexgen)
    orient(a, b, c) ∈ (-1, 0, 1)
end

# Note: This should be true for all, if not this is a bug in the port
@check function orient_against_c(a=complexgen, b=complexgen, c=complexgen)
    dir = orient(a, b, c)
    event!("Port", dir)
    cdir = CPredicates.orient(a, b, c)
    event!("Original", cdir)
    dir == cdir
end

# The next two should agree
check_range(x) = iszero(x) || exponent(abs(x)) ∈ -142:201   # Shewchuk's predicates are only guaranteed to be accurate for floats with binary 
                                                            # exponents in the range [-142, 201]. See Section 2.1 in the paper 
                                                            # Richard Shewchuk, J. Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates. Discrete Comput Geom 18(3), 305–363 (1997)
@check function orient_against_exact(a=complexgen2, b=complexgen2, c=complexgen2)
    assume!(all(check_range, (a.re, a.im, b.re, b.im, c.re, c.im)))
    dir = orient(a, b, c)
    event!("Adaptive", dir)
    edir = ExactPredicates.orient(a, b, c)
    event!("Exact", edir)
    dir == edir
end

@check function orient_c_against_exact(a=complexgen2, b=complexgen2, c=complexgen2)
    assume!(all(check_range, (a.re, a.im, b.re, b.im, c.re, c.im)))
    dir = CPredicates.orient(a, b, c)
    event!("Original", dir)
    edir = ExactPredicates.orient(a, b, c)
    event!("Exact", edir)
    dir == edir
end