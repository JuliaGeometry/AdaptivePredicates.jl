module CPredicates
    const libpredicates = joinpath(@__DIR__, "libpredicates.so")

    function __init__()
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