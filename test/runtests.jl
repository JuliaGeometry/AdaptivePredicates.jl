using AdaptivePredicates
using Test
using Supposition
using BenchmarkTools
import ExactPredicates: ExactPredicates
using Aqua

@testset "Aqua" begin
    Aqua.test_all(AdaptivePredicates)
end

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

setup_orient2(T) = ntuple(_ -> (_rand(T), _rand(T)), 3)
setup_orient3(T) = ntuple(_ -> (_rand(T), _rand(T), _rand(T)), 4)
setup_incircle(T) = ntuple(_ -> (_rand(T), _rand(T)), 4)
setup_insphere(T) = ntuple(_ -> (_rand(T), _rand(T), _rand(T)), 5)
@testset "Allocations" begin
    @test iszero(@ballocated orient2(args...) setup = (args = setup_orient2(Float64)))
    @test iszero(@ballocated orient2(args...) setup = (args = setup_orient2(Float32)))
    @test iszero(@ballocated orient3(args...) setup = (args = setup_orient3(Float64)))
    @test iszero(@ballocated orient3(args...) setup = (args = setup_orient3(Float32)))
    @test iszero(@ballocated incircle(args...) setup = (args = setup_incircle(Float64)))
    @test iszero(@ballocated incircle(args...) setup = (args = setup_incircle(Float32)))
    @test iszero(@ballocated insphere(args...) setup = (args = setup_insphere(Float64)))
    @test iszero(@ballocated insphere(args...) setup = (args = setup_insphere(Float32)))
end

## caches 
function check_non_overlapping(args...)
    ranges = first.(parentindices.(args))
    flag = !issorted(ranges)
    flag && return false
    for i in 1:(length(ranges)-1)
        rᵢ, rᵢ₊₁ = ranges[i], ranges[i+1]
        !isdisjoint(rᵢ, rᵢ₊₁) && return false
    end
    return true
end
function check_gaps(args...)
    ranges = first.(parentindices.(args))
    firstdiffs = diff(collect(first.(ranges)))
    boundaries = firstdiffs .% 64
    return all(iszero, boundaries)
end
v = zeros(100)
v1, v2, v3 = view(v, 1:3), view(v, 2:4), view(v, 7:100)
@test !check_non_overlapping(v1, v2, v3)
@test !check_gaps(v1, v2, v3)

function test_cache(f, lengths) 
    for T in (Float64, Float32)
        args = f(T)
        flag = check_non_overlapping(args...)
        flag = flag && all(lengths .== length.(args))
        flag = flag && all(Base.Fix1(===, parent(args[1])), parent.(args))
        flag = flag && all(x -> eltype(x) === T, args)
        flag = flag && check_gaps(args...)
        return flag
    end
end
@test test_cache(AP.incircleexact_cache, (48, 48, 96, 96, 96, 96, 192, 192, 384))
@test test_cache(AP.incircleslow_cache, (64, 64, 64, 64, 64, 128, 128, 192, 192, 384, 384, 384, 768, 1152))
@test test_cache(AP.incircleadapt_cache, (48, 64, 1152, 1152))





args = AP.incircleexact_cache(Float64)
@test check_non_overlapping(args...)
lengths = (48, 48, 96, 96, 96, 96, 192, 192, 384)
@test all(lengths .== length.(args))
@test all(Base.Fix1(===, parent(args[1])), parent.(args))


regex = r"cache\.([a-zA-Z_][a-zA-Z0-9_]*)"
m = [m.captures[1] for m in eachmatch(regex, str)] |> sort |> unique



using DelaunayTriangulation, BenchmarkTools, Random, Distributions 
const DT = DelaunayTriangulation
struct TriangleSampler{T}
    p::NTuple{2, T}
    q::NTuple{2, T}
    r::NTuple{2, T}
end
Random.eltype(::Type{TriangleSampler{T}}) where {T} = NTuple{2, T}
function Random.rand(rng::Random.AbstractRNG, v::Random.SamplerTrivial{<:TriangleSampler{T}}) where {T}
    # https://blogs.sas.com/content/iml/2020/10/19/random-points-in-triangle.html
    itr = v[]
    p, q, r = itr.p, itr.q, itr.r
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    ax, ay = qx - px, qy - py
    bx, by = rx - px, ry - py
    u₁ = rand(rng, T)
    u₂ = rand(rng, T)
    if u₁ + u₂ > 1
        u₁ = one(T) - u₁
        u₂ = one(T) - u₂
    end
    wx = u₁ * ax + u₂ * bx
    wy = u₁ * ay + u₂ * by
    return (px + wx, py + wy)
end
Random.eltype(::Type{T}) where {P, T<:Triangulation{P}} = NTuple{2, DT.number_type(P)}
function Random.Sampler(::Type{<:Random.AbstractRNG}, tri::Triangulation, ::Random.Repetition)
    V = DT.triangle_type(tri)
    F = DT.number_type(tri)
    triangles = Vector{V}()
    areas = Vector{F}()
    sizehint!(triangles, DT.num_solid_triangles(tri))
    sizehint!(areas, DT.num_solid_triangles(tri))
    for T in DT.each_solid_triangle(tri)
        push!(triangles, T)
        i, j, k = triangle_vertices(T)
        p, q, r = get_point(tri, i, j, k)
        push!(areas, DT.triangle_area(p, q, r))
    end
    A = sum(areas)
    areas ./= A
    return Random.SamplerSimple(tri, (; tri, triangles, areas))
end
function Random.rand(rng::Random.AbstractRNG, spl::Random.SamplerSimple{<:Triangulation})
    tri = spl.data.tri
    triangles = spl.data.triangles
    areas = spl.data.areas
    multi = Multinomial(1, areas)
    T = rand(rng, multi)
    idx = findfirst(isone, T)
    i, j, k = triangle_vertices(triangles[idx])
    p, q, r = get_point(tri, i, j, k)
    return rand(rng, TriangleSampler(p, q, r))
end
function Random.rand!(rng::Random.AbstractRNG, a::AbstractVector, spl::Random.SamplerSimple{<:Triangulation})
    n = length(a)
    tri = spl.data.tri
    triangles = spl.data.triangles
    areas = spl.data.areas
    multi = Multinomial(1, areas)
    samples = rand(rng, multi, n)
    ntri = sum(samples; dims = 2)
    cur_idx = 0
    for (idx, T) in enumerate(triangles)
        ntri[idx] == 0 && continue
        i, j, k = triangle_vertices(T)
        p, q, r = get_point(tri, i, j, k)
        _spl = TriangleSampler(p, q, r)
        m = ntri[idx]
        rand!(rng, view(a, (cur_idx+1):(cur_idx+m)), _spl)
        cur_idx += m
    end
    return a
end
function clear_tls!()
    empty!(task_local_storage())
end
function randpts(n)
    return tuple.(rand(n), rand(n))
end
function collinearpts(nx, ny)
    return DelaunayTriangulation.get_lattice_points(0.0, 1.0, 0.0, 1.0, nx, ny, LinearIndices((1:nx, 1:ny)))
end
function point_location(tri, points)
    return map(points) do p 
        find_triangle(tri, p)
    end
end

function trials()
    n = 1000
    times = zeros(n)
    for i in 1:n 
        @show i
        yield()
        Random.seed!(0)
        clear_tls!()
        points = collinearpts(50, 50)
        times[i] = @elapsed triangulate(points, randomise=false)
    end
    return sum(times) / n
end
ts = trials() # 0.1686563341

#=
BenchmarkTools.Trial: 79 samples with 1 evaluation.
 Range (min … max):  32.860 ms … 120.043 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     60.808 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   63.115 ms ±  21.673 ms  ┊ GC (mean ± σ):  1.90% ± 3.98%

    █▂    ▂ ▂  ▅▂▂    ▂ ▅  ▂    ▂    ▂
  █▁███▁███▅█▁█████▅▅▁████▁█▅█▁▅█▅▅█▅█▁▁▅▅▁▁▁▁▁▁▁█▅▁▁▁▁█▁▅▁▁▁█ ▁
  32.9 ms         Histogram: frequency by time          116 ms <

 Memory estimate: 8.52 MiB, allocs estimate: 23079.
=#

elen, e, flen, f, h = (18, (1.5001756191031835e-22, 6.735097776797844e-16, -1.0497840818857185e-8, -5.960464477539063e-8, 4.603989124298096, -16.0, -8.334350144e9, -4.191888080896e12, -5.670313405836165e17, -1.152921504606847e18, 1.0853497469963175e29, -5.869222279056702e33, 8.582339753016795e45, -4.567192616659072e46, 8.622859660252327e49, 8.62689044908283e62, -5.7066668949226236e79, -5.883942385316154e95, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), 16, (-0.006894129499377399, 0.03125, -1.1030989035655112e14, -8.44424930131968e15, 2.2193738963681804e19, -6.14549659238308e31, 6.003593242680894e34, 3.5434384971605455e47, 2.790737376483359e51, 1.1972621413014757e52, 3.831238852164722e53, 2.2942875581393308e64, -7.351851081381977e68, 1.938453298228514e85, 5.680466505642381e101, 1.0510219116784817e118, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), zeros(32))

