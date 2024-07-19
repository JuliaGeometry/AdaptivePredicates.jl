include("orient2.jl")
include("orient3.jl")
include("incircle.jl")
include("insphere.jl")

_tup(c::Complex) = (c.re, c.im)
sgn(x) = Int(sign(x))

orient2(pa::Complex, pb::Complex, pc::Complex) = orient2(_tup(pa), _tup(pb), _tup(pc))
incircle(pa::Complex, pb::Complex, pc::Complex, pd::Complex) = incircle(_tup(pa), _tup(pb), _tup(pc), _tup(pd))

orient2p(pa, pb, pc) = sgn(orient2(pa, pb, pc))
orient3p(pa, pb, pc, pd) = sgn(orient3(pa, pb, pc, pd))
incirclep(pa, pb, pc, pd) = sgn(incircle(pa, pb, pc, pd))
inspherep(pa, pb, pc, pd, pe) = sgn(insphere(pa, pb, pc, pd, pe))

@doc """
    orient2(p, q, r)
    
For three points `p`, `q`, `r` in the plane, returns:

- A positive value if `r` is left of the oriented line defined by `(p, q)`.
- A negative value if `r` is right of the line.
- Zero if the three points are collinear.

Equivalently, the returned value is positive if the triangle `(p, q, r)` is oriented 
counter-clockwise; negative if the triangle is oriented clockwise; and zero if the points 
are collinear.

The inputs can be either all `NTuple`s or all `Complex`, with either all `Float64` or all `Float32` coordinates.

Note that `orient2(p, q, r)/2` is also roughly the signed area of the triangle defined by `(p, q, r)`. 

See also [`orient2p`](@ref).
"""
orient2

@doc """
    incircle(p, q, r, s)

For four points `p`, `q`, `r`, `s` in the plane, returns:

- A positive value if `s` is inside the circle through `p`, `q`, and `r`.
- A negative value if `s` is outside the circle.
- Zero if the four points are cocircular.

The points `p`, `q`, and `r` are assumed to be in counter-clockwise order, otherwise 
the sign of the result will be reversed.

See also [`incirclep`](@ref).
"""
incircle 

@doc """
    incirclep(p, q, r, s) -> Int 

This function is the same as [`incircle`], except returns the sign of the result,
i.e. `Int(sign(incircle(p, q, r, s)))`.
"""
incirclep

@doc """
    insphere(p, q, r, s, t)

For five points `p`, `q`, `r`, `s`, `t` in the plane, returns: 

- A positive value if `t` is inside the sphere through `p`, `q`, `r`, and `s`.
- A negative value if `t` is outside the sphere.
- Zero if the five points are cospherical. 

The points `p`, `q`, `r`, and `s` are assumed to be ordered so that they 
have positive oriented (according to [`orient3`](@ref)), otherwise 
the sign of the result will be reversed.

See also [`inspherep`](@ref).
"""
insphere

@doc """
    inspherep(p, q, r, s, t) -> Int 

This function is the same as [`insphere`], except returns the sign of the result,
i.e. `Int(sign(incircle(p, q, r, s, t)))`.
"""
inspherep

@doc """
    orient2p(p, q, r) -> Int
    
This function is the same as [`orient2`], except returns the sign of the result,
i.e. `Int(sign(orient2(p, q, r)))`.
"""
orient2p

@doc """
    orient3(p, q, r, s)
    
For four points `p`, `q`, `r`, `s` in space, consider the oriented plane 
on which the triangle `(p, q, r)` is positively oriented. The function returns 

- A positive value if `s` is below this plane.
- A negative value if `s` is below this plane.
- Zero if the four points are coplanar.

The inputs should all be `NTuple`s with either all `Float64` or all `Float32` coordinates.

Note that `orient3(p, q, r, s)/2` is also roughly the signed area of the tetrahedron defined by `(p, q, r, s)`. 

The result returned from this predicate can also be interpreted in terms of the orientation of the tetrahedron 
defined by `(p, q, r, s)`, where a positively oriented tetrahedron `(p, q, r, s)` takes the form

                                   z.
                                 .
                               ,/
                             s
                           ,/|'\\
                         ,/  |  '\\
                       ,/    '.   '\\
                     ,/       |     '\\                 
                   ,/         |       '\\              
                  p-----------'.--------q --> x
                   '\\.         |      ,/              
                      '\\.      |    ,/                 
                         '\\.   '. ,/    
                            '\\. |/      
                               'r       
                                 '\\.    
"""
orient3

@doc """
    orient3p(p, q, r, s) -> Int
    
This function is the same as [`orient3`], except returns the sign of the result,
i.e. `Int(sign(orient2(p, q, r, s)))`.
"""
orient3p