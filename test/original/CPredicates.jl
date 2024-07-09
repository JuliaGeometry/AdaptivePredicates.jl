module CPredicates
const libpredicates = joinpath(@__DIR__, "libpredicates.so")

function _init_()
    exactinit()
end

function orient2d(a, b, c)
    @ccall libpredicates.orient2d(
        Ref(a)::Ptr{Float64}, Ref(b)::Ptr{Float64}, Ref(c)::Ptr{Float64})::Float64
end
function orient2dfast(a, b, c)
    @ccall libpredicates.orient2dfast(
        Ref(a)::Ptr{Float64}, Ref(b)::Ptr{Float64}, Ref(c)::Ptr{Float64})::Float64
end
function orient2dexact(a, b, c)
    @ccall libpredicates.orient2dexact(
        Ref(a)::Ptr{Float64}, Ref(b)::Ptr{Float64}, Ref(c)::Ptr{Float64})::Float64
end
function orient2adapt(a, b, c, detsum)
    @ccall libpredicates.orient2dexact(
        Ref(a)::Ptr{Float64}, Ref(b)::Ptr{Float64}, Ref(c)::Ptr{Float64}, detsum::Float64)::Float64
end
function exactinit()
    @ccall libpredicates.exactinit()::Cvoid
end

function orient(a::Complex, b::Complex, c::Complex)
    orient((a.re, a.im), (b.re, b.im), (c.re, c.im))
end

function orient(a, b, c)
    dir = orient2d(a, b, c)
    if dir < 0.0
        return -1
    elseif dir > 0.0
        return 1
    else
        return 0
    end
end

function orient_exact(a::Complex, b::Complex, c::Complex)
    orient_exact((a.re, a.im), (b.re, b.im), (c.re, c.im))
end

function orient_exact(a, b, c)
    dir = orient2dexact(a, b, c)
    if dir < 0.0
        return -1
    elseif dir > 0.0
        return 1
    else
        return 0
    end
end
end

using .CPredicates: libpredicates
_epsilon,
_splitter,
_resulterrbound,
_ccwerrboundA,
_ccwerrboundB,
_ccwerrboundC,
_o3derrboundA,
_o3derrboundB,
_o3derrboundC,
_iccerrboundA,
_iccerrboundB,
_iccerrboundC,
_isperrboundA,
_isperrboundB,
_isperrboundC = @ccall libpredicates.exactinit2()::NTuple{15,Cdouble}
const AP = AdaptivePredicates
@test _epsilon == AP.epsilon
@test _splitter == AP.splitter
@test _resulterrbound == AP.resulterrbound
@test _ccwerrboundA == AP.ccwerrboundA
@test _ccwerrboundB == AP.ccwerrboundB
@test _ccwerrboundC == AP.ccwerrboundC
@test _o3derrboundA == AP.o3derrboundA
@test _o3derrboundB == AP.o3derrboundB
@test _o3derrboundC == AP.o3derrboundC
@test _iccerrboundA == AP.iccerrboundA
@test _iccerrboundB == AP.iccerrboundB
@test _iccerrboundC == AP.iccerrboundC
@test _isperrboundA == AP.isperrboundA
@test _isperrboundB == AP.isperrboundB
@test _isperrboundC == AP.isperrboundC

const F64 = Float64
const NTESTS = 1000
const L = 1.793662034335766e-43
const R = 3.2138760885179806e10
function _rand()
    return rand((-1, 1)) * rand() * 2.0^(rand(-50:50)) # hard to get high quality floats with just rand() that cover a sufficent distribution, this will do
end
function _rand(n)
    return [_rand() for _ in 1:n]
end

@testset "Test macros" begin
    for _ in 1:NTESTS
        a = _rand()
        @test AdaptivePredicates.Absolute(a) == @ccall libpredicates._Absolute(a::F64)::F64

        a, b, x = _rand(3)
        @test AdaptivePredicates.Fast_Two_Sum_Tail(a, b, x) ==
              @ccall libpredicates._Fast_Two_Sum_Tail(a::F64, b::F64, x::F64)::F64

        a, b = _rand(2)
        @test AdaptivePredicates.Fast_Two_Sum(a, b) ==
              @ccall libpredicates._Fast_Two_Sum(a::F64, b::F64)::NTuple{2,Cdouble}

        a, b, x = _rand(3)
        @test AdaptivePredicates.Fast_Two_Diff_Tail(a, b, x) ==
              @ccall libpredicates._Fast_Two_Diff_Tail(a::F64, b::F64, x::F64)::F64

        a, b = _rand(2)
        @test AdaptivePredicates.Fast_Two_Diff(a, b) ==
              @ccall libpredicates._Fast_Two_Diff(a::F64, b::F64)::NTuple{2,Cdouble}

        a, b, x = _rand(3)
        @test AdaptivePredicates.Two_Sum_Tail(a, b, x) ==
              @ccall libpredicates._Two_Sum_Tail(a::F64, b::F64, x::F64)::F64

        a, b = _rand(2)
        @test AdaptivePredicates.Two_Sum(a, b) ==
              @ccall libpredicates._Two_Sum(a::F64, b::F64)::NTuple{2,Cdouble}

        a, b, x = _rand(3)
        @test AdaptivePredicates.Two_Diff_Tail(a, b, x) ==
              @ccall libpredicates._Two_Diff_Tail(a::F64, b::F64, x::F64)::F64

        a, b = _rand(2)
        @test AdaptivePredicates.Two_Diff(a, b) ==
              @ccall libpredicates._Two_Diff(a::F64, b::F64)::NTuple{2,Cdouble}

        a = _rand()
        @test AdaptivePredicates.Split(a) ==
              @ccall libpredicates._Split(a::F64)::NTuple{2,Cdouble}

        a, b, x = _rand(3)
        @test AdaptivePredicates.Two_Product_Tail(a, b, x) ==
              @ccall libpredicates._Two_Product_Tail(a::F64, b::F64, x::F64)::F64

        a, b = _rand(2)
        @test AdaptivePredicates.Two_Product(a, b) ==
              @ccall libpredicates._Two_Product(a::F64, b::F64)::NTuple{2,Cdouble}

        a, b, bhi, blo = _rand(4)
        @test AdaptivePredicates.Two_Product_Presplit(a, b, bhi, blo) ==
              @ccall libpredicates._Two_Product_Presplit(a::F64, b::F64, bhi::F64, blo::F64)::NTuple{2,Cdouble}

        a, b, ahi, alo, bhi, blo = _rand(6)
        @test AdaptivePredicates.Two_Product_2Presplit(a, ahi, alo, b, bhi, blo) ==
              @ccall libpredicates._Two_Product_2Presplit(a::F64, ahi::F64, alo::F64,
            b::F64, bhi::F64, blo::F64)::NTuple{2,Cdouble}

        a, x = _rand(2)
        @test AdaptivePredicates.Square_Tail(a, x) ==
              @ccall libpredicates._Square_Tail(a::F64, x::F64)::F64

        a = _rand()
        @test AdaptivePredicates.Square(a) ==
              @ccall libpredicates._Square(a::F64)::NTuple{2,Cdouble}

        a1, a0, b = _rand(3)
        @test AdaptivePredicates.Two_One_Sum(a1, a0, b) ==
              @ccall libpredicates._Two_One_Sum(a1::F64, a0::F64, b::F64)::NTuple{3,Cdouble}

        a1, a0, b = _rand(3)
        @test AdaptivePredicates.Two_One_Diff(a1, a0, b) ==
              @ccall libpredicates._Two_One_Diff(a1::F64, a0::F64, b::F64)::NTuple{3,Cdouble}

        a1, a0, b1, b0 = _rand(4)
        @test AdaptivePredicates.Two_Two_Sum(a1, a0, b1, b0) ==
              @ccall libpredicates._Two_Two_Sum(a1::F64, a0::F64, b1::F64, b0::F64)::NTuple{4,Cdouble}

        a1, a0, b1, b0 = _rand(4)
        @test AdaptivePredicates.Two_Two_Diff(a1, a0, b1, b0) ==
              @ccall libpredicates._Two_Two_Diff(a1::F64, a0::F64, b1::F64, b0::F64)::NTuple{4,Cdouble}

        a3, a2, a1, a0, b = _rand(5)
        @test AdaptivePredicates.Four_One_Sum(a3, a2, a1, a0, b) ==
              @ccall libpredicates._Four_One_Sum(a3::F64, a2::F64, a1::F64, a0::F64, b::F64)::NTuple{5,Cdouble}

        a3, a2, a1, a0, b1, b0 = _rand(6)
        @test AdaptivePredicates.Four_Two_Sum(a3, a2, a1, a0, b1, b0) ==
              @ccall libpredicates._Four_Two_Sum(a3::F64, a2::F64, a1::F64, a0::F64, b1::F64, b0::F64)::NTuple{6,Cdouble}

        a3, a2, a1, a0, b4, b3, b1, b0 = _rand(8)
        @test AdaptivePredicates.Four_Four_Sum(a3, a2, a1, a0, b4, b3, b1, b0) ==
              @ccall libpredicates._Four_Four_Sum(a3::F64, a2::F64, a1::F64, a0::F64, b4::F64, b3::F64, b1::F64, b0::F64)::NTuple{8,Cdouble}

        a7, a6, a5, a4, a3, a2, a1, a0, b = _rand(9)
        @test AdaptivePredicates.Eight_One_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b) ==
              @ccall libpredicates._Eight_One_Sum(a7::F64, a6::F64, a5::F64, a4::F64, a3::F64, a2::F64, a1::F64, a0::F64, b::F64)::NTuple{9,Cdouble}

        a7, a6, a5, a4, a3, a2, a1, a0, b1, b0 = _rand(10)
        @test AdaptivePredicates.Eight_Two_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b1, b0) ==
              @ccall libpredicates._Eight_Two_Sum(a7::F64, a6::F64, a5::F64, a4::F64, a3::F64, a2::F64, a1::F64, a0::F64, b1::F64, b0::F64)::NTuple{10,Cdouble}

        a7, a6, a5, a4, a3, a2, a1, a0, b4, b3, b1, b0 = _rand(12)
        @test AdaptivePredicates.Eight_Four_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b4, b3, b1, b0) ==
              @ccall libpredicates._Eight_Four_Sum(a7::F64, a6::F64, a5::F64, a4::F64, a3::F64, a2::F64, a1::F64,
            a0::F64, b4::F64, b3::F64, b1::F64, b0::F64)::NTuple{12,Cdouble}

        a1, a0, b = _rand(3)
        @test AdaptivePredicates.Two_One_Product(a1, a0, b) ==
              @ccall libpredicates._Two_One_Product(a1::F64, a0::F64, b::F64)::NTuple{4,Cdouble}

        a3, a2, a1, a0, b = _rand(5)
        @test AdaptivePredicates.Four_One_Product(a3, a2, a1, a0, b) ==
              @ccall libpredicates._Four_One_Product(a3::F64, a2::F64, a1::F64, a0::F64, b::F64)::NTuple{8,Cdouble}

        a3, a2, a1, a0 = _rand(4)
        @test AdaptivePredicates.Two_Two_Product(a3, a2, a1, a0) ==
              @ccall libpredicates._Two_Two_Product(a3::F64, a2::F64, a1::F64, a0::F64)::NTuple{8,Cdouble}

        a1, a0 = _rand(2)
        @test collect(AdaptivePredicates.Two_Square(a1, a0)) ≈ # Why isn't Two_Square(a1, a0) exactly the same???
              collect(@ccall libpredicates._Two_Square(a1::F64, a0::F64)::NTuple{6,Cdouble}) atol = 1e-14 rtol = 1e-10
    end
end