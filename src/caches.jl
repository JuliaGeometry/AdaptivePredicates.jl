@static if isdefined(Base, :Memory)
    const Vec{T} = Memory{T}
else
    const Vec{T} = Vector{T}
end

struct CacheKey{T}
    size::UInt16
    id::UInt8
end

struct APMarker{T} end # Avoid collisions with other packages using the task_local_storage

const APCache{T} = Dict{CacheKey{T},Vec{T}}

@inline function task_local_cache(::Type{T}) where {T}
    tls = get!(task_local_storage(), APMarker{T}()) do
        APCache{T}()
    end::APCache{T}
    return tls::APCache{T}
end


@inline function get_cache!(tls::APCache{T}, size, id) where {T}
    cache::Vec{T} = get!(tls, CacheKey{T}(size, id)) do
        Vec{T}(undef, Int(size)) # Memory{T}(undef, Int(size)) has weird concurrency issues sometimes?
    end
    return cache::Vec{T}
end

struct Orient3Cache{T}
    h4::NTuple{4,T}
    h8::NTuple{8,T}
    h12::NTuple{12,T}
    h16::NTuple{16,T}
    h24::NTuple{24,T}
    h32::NTuple{32,T}
    h48_1::Vec{T}
    h48_2::Vec{T}
    h64_1::Vec{T}
    h64_2::Vec{T}
    h64_3::Vec{T}
    h96::Vec{T}
    h128::Vec{T}
    h192_1::Vec{T}
    h192_2::Vec{T}
    h196::Vec{T}
end
@inline function Orient3Cache{T}() where {T}
    tls = task_local_cache(T)
    return Orient3Cache{T}(tls)
end
@inline function Orient3Cache{T}(tls::APCache{T}) where {T}
    h4 = ntuple(_ -> zero(T), Val(4))
    h8 = ntuple(_ -> zero(T), Val(8))
    h12 = ntuple(_ -> zero(T), Val(12))
    h16 = ntuple(_ -> zero(T), Val(16))
    h24 = ntuple(_ -> zero(T), Val(24))
    h32 = ntuple(_ -> zero(T), Val(32))
    h48_1 = get_cache!(tls, 0x0030, 0x01)::Vec{T}
    h48_2 = get_cache!(tls, 0x0030, 0x02)::Vec{T}
    h64_1 = get_cache!(tls, 0x0040, 0x01)::Vec{T}
    h64_2 = get_cache!(tls, 0x0040, 0x02)::Vec{T}
    h64_3 = get_cache!(tls, 0x0040, 0x03)::Vec{T}
    h96 = get_cache!(tls, 0x0060, 0x01)::Vec{T}
    h128 = get_cache!(tls, 0x0080, 0x01)::Vec{T}
    h192_1 = get_cache!(tls, 0x00c0, 0x01)::Vec{T}
    h192_2 = get_cache!(tls, 0x00c0, 0x02)::Vec{T}
    h196 = get_cache!(tls, 0x00c4, 0x01)::Vec{T}
    return Orient3Cache{T}(
        h4, h8, h12, h16, h24, h32, h48_1, h48_2, h64_1, h64_2, h64_3,
        h96, h128, h192_1, h192_2, h196
    )
end

function incircleexact_cache(::Type{T}) where {T}
    cache_size = 1408                   # (64 + 64) + (128 + 128 + 128 + 128) + (192 + 192) + (384)
    cache = Vec{T}(undef, cache_size)
    h48_1 = view(cache, 1:48)           # 0 .+ (1:48)
    h48_2 = view(cache, 65:112)         # 64 .+ (1:48)
    h96_1 = view(cache, 129:224)        # 128 .+ (1:96)
    h96_2 = view(cache, 257:352)        # 256 .+ (1:96)
    h96_3 = view(cache, 385:480)        # 384 .+ (1:96)
    h96_4 = view(cache, 513:608)        # 512 .+ (1:96)
    h192_1 = view(cache, 641:832)       # 640 .+ (1:192)
    h192_2 = view(cache, 833:1024)      # 832 .+ (1:192)
    h384_1 = view(cache, 1025:1408)     # 1024 .+ (1:384)
    return h48_1, h48_2, h96_1, h96_2, h96_3, h96_4, h192_1, h192_2, h384_1
end
function incircleslow_cache(::Type{T}) where {T}
    cache_size = (64 + 64 + 64 + 64 + 64 + 64) + (128 + 128) + (192 + 192) + (384 + 384 + 384) + (768) + (1152)
    cache = Vec{T}(undef, cache_size)
    offset = 0
    h64_1 = view(cache, 0 .+ (1:64))
    offset += 64
    h64_2 = view(cache, offset .+ (1:64))
    offset += 64
    h64_3 = view(cache, offset .+ (1:64))
    offset += 64
    h64_4 = view(cache, offset .+ (1:64))
    offset += 64
    h64_5 = view(cache, offset .+ (1:64))
    offset += 64
    h64_6 = view(cache, offset .+ (1:64))
    offset += 64
    h128_1 = view(cache, offset .+ (1:128))
    offset += 128
    h128_2 = view(cache, offset .+ (1:128))
    offset += 128
    h192_1 = view(cache, offset .+ (1:192))
    offset += 192
    h192_2 = view(cache, offset .+ (1:192))
    offset += 192
    h384_1 = view(cache, offset .+ (1:384))
    offset += 384
    h384_2 = view(cache, offset .+ (1:384))
    offset += 384
    h384_3 = view(cache, offset .+ (1:384))
    offset += 384
    h768 = view(cache, offset .+ (1:768))
    offset += 768
    h1152_1 = view(cache, offset .+ (1:1152))
    return h64_1, h64_2, h64_3, h64_4, h64_5, h64_6, h128_1, h128_2, h192_1, h192_2, h384_1, h384_2, h384_3, h768, h1152_1
end
incircleadapt_cache(_, cache) = cache
incircleadapt_cache(::Type{T}, ::Nothing) where {T} = incircleadapt_cache(T)
function incircleadapt_cache(::Type{T}) where {T}
    @inbounds begin
        cache_size = (64) + (64) + (1152 + 1152)
        cache = Vec{T}(undef, cache_size)
        offset = 0
        h48_1 = view(cache, 0 .+ (1:48))
        offset += 64
        h64_1 = view(cache, offset .+ (1:64))
        offset += 64
        h1152_1 = view(cache, offset .+ (1:1152))
        offset += 1152
        h1152_2 = view(cache, offset .+ (1:1152))
        return h48_1, h64_1, h1152_1, h1152_2
    end
end

struct InsphereCache{T}
    h4::NTuple{4,T}
    h8::NTuple{8,T}
    h12::NTuple{12,T}
    h16::NTuple{16,T}
    h24::NTuple{24,T}
    h32::NTuple{32,T}
    h48_1::Vec{T}
    h48_2::Vec{T}
    h64_1::Vec{T}
    h64_2::Vec{T}
    h64_3::Vec{T}
    h96_1::Vec{T}
    h96_2::Vec{T}
    h96_3::Vec{T}
    h96_4::Vec{T}
    h96_5::Vec{T}
    h128::Vec{T}
    h192::Vec{T}
    h288_1::Vec{T}
    h288_2::Vec{T}
    h288_3::Vec{T}
    h288_4::Vec{T}
    h384_1::Vec{T}
    h384_2::Vec{T}
    h384_3::Vec{T}
    h384_4::Vec{T}
    h384_5::Vec{T}
    h384_6::Vec{T}
    h576_1::Vec{T}
    h576_2::Vec{T}
    h768_1::Vec{T}
    h768_2::Vec{T}
    h768_3::Vec{T}
    h768_4::Vec{T}
    h768_5::Vec{T}
    h768_6::Vec{T}
    h768_7::Vec{T}
    h768_8::Vec{T}
    h768_9::Vec{T}
    h1152_1::Vec{T}
    h1152_2::Vec{T}
    h1152_3::Vec{T}
    h1152_4::Vec{T}
    h1152_5::Vec{T}
    h1536_1::Vec{T}
    h1536_2::Vec{T}
    h1536_3::Vec{T}
    h2304_1::Vec{T}
    h2304_2::Vec{T}
    h2304_3::Vec{T}
    h3456::Vec{T}
    h4608::Vec{T}
    h5760::Vec{T}
    h6912_1::Vec{T}
    h6912_2::Vec{T}
    h6912_3::Vec{T}
    h6912_4::Vec{T}
    h13824_1::Vec{T}
    h13824_2::Vec{T}
    h27648::Vec{T}
end
@inline function InsphereCache{T}() where {T}
    tls = task_local_cache(T)
    return InsphereCache{T}(tls)
end
@inline function InsphereCache{T}(tls::APCache{T}) where {T}
    h4 = ntuple(_ -> zero(T), Val(4))
    h8 = ntuple(_ -> zero(T), Val(8))
    h12 = ntuple(_ -> zero(T), Val(12))
    h16 = ntuple(_ -> zero(T), Val(16))
    h24 = ntuple(_ -> zero(T), Val(24))
    h32 = ntuple(_ -> zero(T), Val(32))
    h48_1 = get_cache!(tls, 0x0030, 0x01)::Vec{T}
    h48_2 = get_cache!(tls, 0x0030, 0x02)::Vec{T}
    h64_1 = get_cache!(tls, 0x0040, 0x01)::Vec{T}
    h64_2 = get_cache!(tls, 0x0040, 0x02)::Vec{T}
    h64_3 = get_cache!(tls, 0x0040, 0x03)::Vec{T}
    h96_1 = get_cache!(tls, 0x0060, 0x01)::Vec{T}
    h96_2 = get_cache!(tls, 0x0060, 0x02)::Vec{T}
    h96_3 = get_cache!(tls, 0x0060, 0x03)::Vec{T}
    h96_4 = get_cache!(tls, 0x0060, 0x04)::Vec{T}
    h96_5 = get_cache!(tls, 0x0060, 0x05)::Vec{T}
    h128 = get_cache!(tls, 0x0080, 0x01)::Vec{T}
    h192 = get_cache!(tls, 0x00c0, 0x01)::Vec{T}
    h288_1 = get_cache!(tls, 0x0120, 0x01)::Vec{T}
    h288_2 = get_cache!(tls, 0x0120, 0x02)::Vec{T}
    h288_3 = get_cache!(tls, 0x0120, 0x03)::Vec{T}
    h288_4 = get_cache!(tls, 0x0120, 0x04)::Vec{T}
    h384_1 = get_cache!(tls, 0x0180, 0x01)::Vec{T}
    h384_2 = get_cache!(tls, 0x0180, 0x02)::Vec{T}
    h384_3 = get_cache!(tls, 0x0180, 0x03)::Vec{T}
    h384_4 = get_cache!(tls, 0x0180, 0x04)::Vec{T}
    h384_5 = get_cache!(tls, 0x0180, 0x05)::Vec{T}
    h384_6 = get_cache!(tls, 0x0180, 0x06)::Vec{T}
    h576_1 = get_cache!(tls, 0x0240, 0x01)::Vec{T}
    h576_2 = get_cache!(tls, 0x0240, 0x02)::Vec{T}
    h768_1 = get_cache!(tls, 0x0300, 0x01)::Vec{T}
    h768_2 = get_cache!(tls, 0x0300, 0x02)::Vec{T}
    h768_3 = get_cache!(tls, 0x0300, 0x03)::Vec{T}
    h768_4 = get_cache!(tls, 0x0300, 0x04)::Vec{T}
    h768_5 = get_cache!(tls, 0x0300, 0x05)::Vec{T}
    h768_6 = get_cache!(tls, 0x0300, 0x06)::Vec{T}
    h768_7 = get_cache!(tls, 0x0300, 0x07)::Vec{T}
    h768_8 = get_cache!(tls, 0x0300, 0x08)::Vec{T}
    h768_9 = get_cache!(tls, 0x0300, 0x09)::Vec{T}
    h1152_1 = get_cache!(tls, 0x0480, 0x01)::Vec{T}
    h1152_2 = get_cache!(tls, 0x0480, 0x02)::Vec{T}
    h1152_3 = get_cache!(tls, 0x0480, 0x03)::Vec{T}
    h1152_4 = get_cache!(tls, 0x0480, 0x04)::Vec{T}
    h1152_5 = get_cache!(tls, 0x0480, 0x05)::Vec{T}
    h1536_1 = get_cache!(tls, 0x0600, 0x01)::Vec{T}
    h1536_2 = get_cache!(tls, 0x0600, 0x02)::Vec{T}
    h1536_3 = get_cache!(tls, 0x0600, 0x03)::Vec{T}
    h2304_1 = get_cache!(tls, 0x0900, 0x01)::Vec{T}
    h2304_2 = get_cache!(tls, 0x0900, 0x02)::Vec{T}
    h2304_3 = get_cache!(tls, 0x0900, 0x03)::Vec{T}
    h3456 = get_cache!(tls, 0x0d80, 0x01)::Vec{T}
    h4608 = get_cache!(tls, 0x1200, 0x01)::Vec{T}
    h5760 = get_cache!(tls, 0x02400, 0x01)::Vec{T}
    h6912_1 = get_cache!(tls, 0x1b00, 0x01)::Vec{T}
    h6912_2 = get_cache!(tls, 0x1b00, 0x02)::Vec{T}
    h6912_3 = get_cache!(tls, 0x1b00, 0x03)::Vec{T}
    h6912_4 = get_cache!(tls, 0x1b00, 0x04)::Vec{T}
    h13824_1 = get_cache!(tls, 0x3600, 0x01)::Vec{T}
    h13824_2 = get_cache!(tls, 0x3600, 0x02)::Vec{T}
    h27648 = get_cache!(tls, 0x6c00, 0x01)::Vec{T}
    return InsphereCache{T}(
        h4, h8, h12, h16, h24, h32, h48_1, h48_2,
        h64_1, h64_2, h64_3, h96_1, h96_2, h96_3, h96_4, h96_5,
        h128, h192, h288_1, h288_2, h288_3, h288_4,
        h384_1, h384_2, h384_3, h384_4, h384_5, h384_6,
        h576_1, h576_2, h768_1, h768_2, h768_3, h768_4, h768_5, h768_6, h768_7, h768_8, h768_9,
        h1152_1, h1152_2, h1152_3, h1152_4, h1152_5,
        h1536_1, h1536_2, h1536_3, h2304_1, h2304_2, h2304_3,
        h3456, h4608, h5760, h6912_1, h6912_2, h6912_3, h6912_4,
        h13824_1, h13824_2, h27648
    )
end