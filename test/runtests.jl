using AdaptivePredicates
using Test
using Supposition
import ExactPredicates: ExactPredicates
const AP = AdaptivePredicates
cd("original") do
    include("compile.jl")
end
include("CPredicates.jl")
import .CPredicates as C

include("utils.jl")

@testset "Check exactinit()" begin
    @test AP.InitConstants{Float64}(C.C64_consts...) == AP.IC64
    @test AP.InitConstants{Float32}(C.C32_consts...) == AP.IC32
end

const MACROS = (
    MacroMethod(:Absolute, 1),
    MacroMethod(:Fast_Two_Sum_Tail, 3),
    MacroMethod(:Fast_Two_Sum, 2),
    MacroMethod(:Fast_Two_Diff_Tail, 3),
    MacroMethod(:Fast_Two_Diff, 2),
    MacroMethod(:Two_Sum_Tail, 3),
    MacroMethod(:Two_Sum, 2),
    MacroMethod(:Two_Diff_Tail, 3),
    MacroMethod(:Two_Diff, 2),
    MacroMethod(:Split, 1),
    MacroMethod(:Two_Product_Tail, 3),
    MacroMethod(:Two_Product, 2),
    MacroMethod(:Two_Product_Presplit, 4),
    MacroMethod(:Two_Product_2Presplit, 6),
    MacroMethod(:Square_Tail, 2),
    MacroMethod(:Square, 1),
    MacroMethod(:Two_One_Sum, 3),
    MacroMethod(:Two_One_Diff, 3),
    MacroMethod(:Two_Two_Sum, 4),
    MacroMethod(:Two_Two_Diff, 4),
    MacroMethod(:Four_One_Sum, 5),
    MacroMethod(:Four_Two_Sum, 6),
    MacroMethod(:Four_Four_Sum, 8),
    MacroMethod(:Eight_One_Sum, 9),
    MacroMethod(:Eight_Two_Sum, 10),
    MacroMethod(:Eight_Four_Sum, 12),
    MacroMethod(:Two_One_Product, 3),
    MacroMethod(:Four_One_Product, 5),
    MacroMethod(:Two_Two_Product, 4),
    MacroMethod(:Two_Square, 2)
)
_sum_size(args) = length(args[end])
_grow_size(args) = length(args[end]) + length(args[end-1])
_scale_size(args) = 2length(args[end-1])
const ARITHMETIC = (
    ArithmeticMethod(:grow_expansion, 3, (-1, 1, _grow_size)),
    ArithmeticMethod(:grow_expansion_zeroelim, 3, (-1, 1, _grow_size)),
    ArithmeticMethod(:expansion_sum, 3, (-1, -1, _grow_size)),
    ArithmeticMethod(:expansion_sum_zeroelim1, 3, (-1, -1, _grow_size)),
    ArithmeticMethod(:expansion_sum_zeroelim2, 3, (-1, -1, _grow_size)),
    ArithmeticMethod(:fast_expansion_sum, 3, (-1, -1, _grow_size)),
    ArithmeticMethod(:fast_expansion_sum_zeroelim, 3, (-1, -1, _grow_size)),
    ArithmeticMethod(:linear_expansion_sum, 3, (-1, -1, _grow_size)),
    ArithmeticMethod(:linear_expansion_sum_zeroelim, 3, (-1, -1, _grow_size)),
    ArithmeticMethod(:scale_expansion, 3, (-1, 1, _scale_size)),
    ArithmeticMethod(:scale_expansion_zeroelim, 3, (-1, 1, _scale_size)),
    ArithmeticMethod(:compress, 2, (-1, _sum_size)),
    ArithmeticMethod(:estimate, 1, (-1,))
)
const PREDICATES = (
    PredicateMethod(:orient2, 3, 2),
    PredicateMethod(:orient3, 4, 3),
    PredicateMethod(:incircle, 4, 2),
    PredicateMethod(:insphere, 5, 3)
)

@testset "Exact Equality Tests" begin
    @testset "Macros" begin
        foreach(MACROS) do f
            @repeat test_f(f)
        end
    end

    @testset "Arithmetic" begin
        foreach(ARITHMETIC) do f
            @repeat test_f(f)
        end
    end

    @testset "Predicates" begin
        foreach(PREDICATES) do f
            @repeat test_f(f)
        end
    end
end

check_range(x::Float64) = iszero(x) || exponent(abs(x)) ∈ -142:201 # Ranges are explained in the README
check_range(x::Float32) = iszero(x) || exponent(abs(x)) ∈ -24:24 # Don't know the exact range for Float32, it might be [-3, 24] but ExactPredicates doesn't use Float32 anyway for me to confirm.

@testset "Supposition tests" begin
    for (T, Tn) in ((Float64, 64), (Float32, 32))
        fgen = Data.Floats{T}(infs=false, nans=false)
        cgen = @composed _complex(a=fgen, b=fgen) = a + b * im
        r2gen = @composed _tuple(a=fgen, b=fgen) = (a, b)
        r3gen = @composed _tuple(a=fgen, b=fgen, c=fgen) = (a, b, c)
        # Against C
        @check function _orient2(p=cgen, q=cgen, r=cgen)
            assume!(all(check_range, (p.re, p.im, q.re, q.im, r.re, r.im)))
            ap = orient2(p, q, r)
            c = C.orient2d((p.re, p.im), (q.re, q.im), (r.re, r.im))
            event!("AdaptiveOrient2$(Tn)", ap)
            event!("COrient2$(Tn)", c)
            ap == c
        end

        @check function _orient3(p=r3gen, q=r3gen, r=r3gen, s=r3gen)
            assume!(all(check_range, (p..., q..., r..., s...)))
            ap = orient3(p, q, r, s)
            c = C.orient3d(p, q, r, s)
            event!("AdaptiveOrient3$(Tn)", ap)
            event!("COrient3$(Tn)", c)
            ap == c
        end

        @check function _incircle(p=r2gen, q=r2gen, r=r2gen, s=r2gen)
            assume!(all(check_range, (p..., q..., r..., s...)))
            ap = incircle(p, q, r, s)
            c = C.incircle(p, q, r, s)
            event!("AdaptiveIncircle$(Tn)", ap)
            event!("CIncircle$(Tn)", c)
            ap == c
        end

        @check function _insphere(p=r3gen, q=r3gen, r=r3gen, s=r3gen, t=r3gen)
            assume!(all(check_range, (p..., q..., r..., s..., t...)))
            ap = insphere(p, q, r, s, t)
            c = C.insphere(p, q, r, s, t)
            event!("AdaptiveInsphere$(Tn)", ap)
            event!("CInsphere$(Tn)", c)
            ap == c
        end

        # Against ExactPredicates
        if T ≠ Float32
            @check function _orient2p(p=cgen, q=cgen, r=cgen)
                assume!(all(check_range, (p.re, p.im, q.re, q.im, r.re, r.im)))
                ap = orient2p(p, q, r)
                c = ExactPredicates.orient((p.re, p.im), (q.re, q.im), (r.re, r.im))
                event!("AdaptiveOrient2$(Tn)", ap)
                event!("ExactPredicatesOrient2$(Tn)", c)
                ap == c
            end

            @check function _orient3p(p=r3gen, q=r3gen, r=r3gen, s=r3gen)
                assume!(all(check_range, (p..., q..., r..., s...)))
                ap = orient3p(p, q, r, s)
                c = ExactPredicates.orient(p, q, r, s)
                event!("AdaptiveOrient3$(Tn)", ap)
                event!("ExactPredicatesOrient3$(Tn)", c)
                ap == c
            end

            @check function _incirclep(p=r2gen, q=r2gen, r=r2gen, s=r2gen)
                assume!(all(check_range, (p..., q..., r..., s...)))
                ap = incirclep(p, q, r, s)
                c = ExactPredicates.incircle(p, q, r, s)
                event!("AdaptiveIncircle$(Tn)", ap)
                event!("ExactPredicatesIncircle$(Tn)", c)
                ap == c
            end

            @check function _inspherep(p=r3gen, q=r3gen, r=r3gen, s=r3gen, t=r3gen)
                assume!(all(check_range, (p..., q..., r..., s..., t...)))
                ap = inspherep(p, q, r, s, t)
                c = ExactPredicates.insphere(p, q, r, s, t)
                event!("AdaptiveInsphere$(Tn)", ap)
                event!("ExactPredicatesInsphere$(Tn)", c)
                ap == c
            end
        end
    end
end

@testset "Allocations" begin
    setup_orient2(T) = ntuple(_ -> (_rand(T), _rand(T)), 3)
    setup_orient3(T) = ntuple(_ -> (_rand(T), _rand(T), _rand(T)), 4)
    setup_incircle(T) = ntuple(_ -> (_rand(T), _rand(T)), 4)
    setup_insphere(T) = ntuple(_ -> (_rand(T), _rand(T), _rand(T)), 5)

    @test iszero(@ballocated orient2(args...) setup = (args = setup_orient2(Float64)))
    @test iszero(@ballocated orient2(args...) setup = (args = setup_orient2(Float32)))
    @test iszero(@ballocated orient3(args...) setup = (args = setup_orient3(Float64)))
    @test iszero(@ballocated orient3(args...) setup = (args = setup_orient3(Float32)))
    @test iszero(@ballocated incircle(args...) setup = (args = setup_incircle(Float64)))
    @test iszero(@ballocated incircle(args...) setup = (args = setup_incircle(Float32)))
    @test iszero(@ballocated insphere(args...) setup = (args = setup_insphere(Float64)))
    @test iszero(@ballocated insphere(args...) setup = (args = setup_insphere(Float32)))
end