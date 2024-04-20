module AdaptivePredicates

using StaticArrays

export orient2

function find_epsilon()
    every_other = true
    half = 0.5
    epsilon = 1.0
    splitter = 1.0
    check = 1.0
    # /* Repeatedly divide `epsilon' by two until it is too small to add to    */
    # /*   one without causing roundoff.  (Also check if the sum is equal to   */
    # /*   the previous sum, for machines that round up instead of using exact */
    # /*   rounding.  Not that this library will work on such machines anyway. */
    cond = true
    while cond
    lastcheck = check
    epsilon *= half
    if (every_other)
        splitter *= 2.0
    end
    every_other = !every_other
    check = 1.0 + epsilon
    cond = ((check != 1.0) && (check != lastcheck))
    end
    splitter += 1.0
    return epsilon, splitter
end

const epsilon, splitter = find_epsilon()
@assert epsilon == eps(1.0) / 2 
const resulterrbound = (3.0 + 8.0 * epsilon) * epsilon
const ccwerrboundA = (3.0 + 16.0 * epsilon) * epsilon
const ccwerrboundB = (2.0 + 12.0 * epsilon) * epsilon
const ccwerrboundC = (9.0 + 64.0 * epsilon) * epsilon * epsilon
const o3derrboundA = (7.0 + 56.0 * epsilon) * epsilon
const o3derrboundB = (3.0 + 28.0 * epsilon) * epsilon
const o3derrboundC = (26.0 + 288.0 * epsilon) * epsilon * epsilon
const iccerrboundA = (10.0 + 96.0 * epsilon) * epsilon
const iccerrboundB = (4.0 + 48.0 * epsilon) * epsilon
const iccerrboundC = (44.0 + 576.0 * epsilon) * epsilon * epsilon
const isperrboundA = (16.0 + 224.0 * epsilon) * epsilon
const isperrboundB = (5.0 + 72.0 * epsilon) * epsilon
const isperrboundC = (71.0 + 1408.0 * epsilon) * epsilon * epsilon


include("utils.jl")

"""
Return a positive value if the points pa, pb, and pc occur
in counterclockwise order a negative value if they occur
in clockwise order and zero if they are collinear.  The
result is also a rough approximation of twice the signed
area of the triangle defined by the three points.                                                             
"""
function orient2(pa, pb, pc)::Float64
  detleft = (pa[1] - pc[1]) * (pb[2] - pc[2])
  detright = (pa[2] - pc[2]) * (pb[1] - pc[1])
  det = detleft - detright

  if (detleft > 0.0)
    if (detright <= 0.0)
      return det
    else
      detsum = detleft + detright
	end
  elseif (detleft < 0.0)
    if (detright >= 0.0) 
      return det
    else 
      detsum = -detleft - detright
	end
  else
    return det
  end

  errbound = ccwerrboundA * detsum
  if (det >= errbound) || (-det >= errbound)
    return det
  end

  return orient2adapt(pa, pb, pc, detsum)
end

function orient2adapt(pa, pb, pc, detsum::Float64)::Float64
  C1 = MVector{8, Float64}(undef)
  C2 = MVector{12, Float64}(undef)
  D  = MVector{16, Float64}(undef)
	
  acx = (pa[1] - pc[1])
  bcx = (pb[1] - pc[1])
  acy = (pa[2] - pc[2])
  bcy = (pb[2] - pc[2])

  detleft, detlefttail = Two_Product(acx, bcy)
  detright, detrighttail = Two_Product(acy, bcx)

  B3, B2, B1, B0 = Two_Two_Diff(detleft, detlefttail, detright, detrighttail)
  B = (B0, B1, B2, B3)

  det = estimate(4, B)
  errbound = ccwerrboundB * detsum
  if ((det >= errbound) || (-det >= errbound))
    return det
  end

  acxtail = Two_Diff_Tail(pa[1], pc[1], acx)
  bcxtail = Two_Diff_Tail(pb[1], pc[1], bcx)
  acytail = Two_Diff_Tail(pa[2], pc[2], acy)
  bcytail = Two_Diff_Tail(pb[2], pc[2], bcy)

  if ((acxtail == 0.0) && (acytail == 0.0)
      && (bcxtail == 0.0) && (bcytail == 0.0))
    return det
  end

  errbound = ccwerrboundC * detsum + resulterrbound * Absolute(det)
  det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail)
  if ((det >= errbound) || (-det >= errbound))
    return det
  end

  s1, s0 = Two_Product(acxtail, bcy)
  t1, t0 = Two_Product(acytail, bcx)
  u3, u2, u1, u0 = Two_Two_Diff(s1, s0, t1, t0)
  u = (u0, u1, u2, u3)
  C1length = fast_expansion_sum_zeroelim(4, B, 4, u, C1)

  s1, s0 = Two_Product(acx, bcytail)
  t1, t0 = Two_Product(acy, bcxtail)
  u3, u2, u1, u0 = Two_Two_Diff(s1, s0, t1, t0)
  u = (u0, u1, u2, u3)
  C2length = fast_expansion_sum_zeroelim(C1length, C1, 4, u, C2)

  s1, s0 = Two_Product(acxtail, bcytail)
  t1, t0 = Two_Product(acytail, bcxtail)
  u3, u2, u1, u0 = Two_Two_Diff(s1, s0, t1, t0)
  u = (u0, u1, u2, u3)
  Dlength = fast_expansion_sum_zeroelim(C2length, C2, 4, u, D)

  return D[Dlength] # originally Dlength - 1 
end


end # module AdaptivePredicates
