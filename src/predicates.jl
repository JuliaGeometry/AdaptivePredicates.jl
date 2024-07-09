
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
  C1, C1length = fast_expansion_sum_zeroelim(4, B, 4, u, Val(8))

  s1, s0 = Two_Product(acx, bcytail)
  t1, t0 = Two_Product(acy, bcxtail)
  u3, u2, u1, u0 = Two_Two_Diff(s1, s0, t1, t0)
  u = (u0, u1, u2, u3)
  C2, C2length = fast_expansion_sum_zeroelim(C1length, C1, 4, u, Val(12))

  s1, s0 = Two_Product(acxtail, bcytail)
  t1, t0 = Two_Product(acytail, bcxtail)
  u3, u2, u1, u0 = Two_Two_Diff(s1, s0, t1, t0)
  u = (u0, u1, u2, u3)
  D, Dlength = fast_expansion_sum_zeroelim(C2length, C2, 4, u, Val(16))  
  return @inbounds D[Dlength] # originally Dlength - 1 
end