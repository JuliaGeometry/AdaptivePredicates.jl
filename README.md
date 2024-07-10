# AdaptivePredicates

A Julia port of "Routines for Arbitrary Precision Floating-point Arithmetic and Fast Robust Geometric Predicates"
by Jonathan Richard Shewchuk. http://www.cs.cmu.edu/~quake/robust.html

Please note that these predicates are only assumed to be correct for floating point numbers with binary exponents in the 
range $[-142, 201]$. In particular:

> This article does not address issues of overflow and underflow, so I allow the exponent to be an integer in the range  $[-\infty, \infty]$. (Fortunately, many applications have inputs whose exponents fall within a circumscribed range. The four predicates implemented for this article will not overflow nor underflow if their inputs have exponents in the range $[-142, 201]$ and IEEE 754 double precision arithmetic is used.)

- Richard Shewchuk, J. Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates. Discrete Comput Geom 18(3), 305–363 (1997)

If you need predicates outside of this range, ExactPredicates.jl might be preferable.

The package is still in development and is not registered. If you want to use the package, you can do
```julia
using Pkg
Pkg.add("https://github.com/JuliaGeometry/AdaptivePredicates.jl")
using AdaptivePredicates
```

## License
The original code is in the public domain and this Julia port is under the MIT License

## Translated primitives
- orient (2D)
