function Split(a)
  c = (splitter * a)
  abig = (c - a)
  ahi = c - abig
  alo = a - ahi
  return ahi, alo
end

function Two_Product_Tail(a, b, x)
  ahi, alo = Split(a)
  bhi, blo = Split(b)
  err1 = x - (ahi * bhi)
  err2 = err1 - (alo * bhi)
  err3 = err2 - (ahi * blo)
  y = (alo * blo) - err3
  return y
end

function Two_Product(a, b)
	x = a * b
	return x, Two_Product_Tail(a, b, x)
end

function Two_Diff_Tail(a, b, x)
  bvirt = (a - x)
  avirt = x + bvirt
  bround = bvirt - b
  around = a - avirt
  y = around + bround
  return y
end

function Two_Diff(a, b)
  x = (a - b)
  return x, Two_Diff_Tail(a, b, x)
end

function Two_Sum_Tail(a, b, x)
  bvirt = (x - a)
  avirt = x - bvirt
  bround = b - bvirt
  around = a - avirt
  y = around + bround
  return y
end

function Two_Sum(a, b) 
  x = (a + b)
  return x, Two_Sum_Tail(a, b, x)
end

function Fast_Two_Sum_Tail(a, b, x)
    bvirt = x - a
    y = b - bvirt
    return y
end

function Fast_Two_Sum(a, b)
    x = (a + b)
    return x, Fast_Two_Sum_Tail(a, b, x)
end

function Two_One_Diff(a1, a0, b)
	i, x0 = Two_Diff(a0, b)
	x2, x1 = Two_Sum(a1, i)
	return x2, x1, x0
end

function Two_Two_Diff(a1, a0, b1, b0)
  _j, _0, x0 = Two_One_Diff(a1, a0, b0)
  x3, x2, x1 = Two_One_Diff(_j, _0, b1)
  return x3, x2, x1, x0
end

Absolute(a) = ((a) >= 0.0 ? (a) : -(a))

function estimate(elen, e)
  Q = e[1]
  for eindex in 2:elen
    Q += e[eindex]
  end
  return Q
end

"""
    Sets h = e + f.  See the long version of my paper for details. 

	If round-to-even is used (as with IEEE 754), maintains the strongly
	nonoverlapping property.  (That is, if e is strongly nonoverlapping, h
	will be also.)  Does NOT maintain the nonoverlapping or nonadjacent
	properties.                  
"""
@inline function fast_expansion_sum_zeroelim(elen::Int, e, flen::Int, f, h)::Int
  @inbounds begin
  enow = e[1]
  fnow = f[1]
  eindex = findex = 1
  if ((fnow > enow) == (fnow > -enow))
    Q = enow
    enow = e[eindex += 1]
  else
    Q = fnow
    fnow = f[findex += 1]
  end
  hindex = 0
  if ((eindex < elen) && (findex < flen)) # still < since pre-increment
    if ((fnow > enow) == (fnow > -enow))
      Qnew, hh = Fast_Two_Sum(enow, Q)
      enow = e[eindex += 1]
    else
      Qnew, hh = Fast_Two_Sum(fnow, Q)
      fnow = f[findex += 1]
	end
    Q = Qnew
    if hh != 0.0
      h[hindex+=1] = hh
	end

    while (eindex < elen) && (findex < flen)  # still < since pre-increment
      if (fnow > enow) == (fnow > -enow) 
        Qnew, hh = Two_Sum(Q, enow)
        enow = e[eindex += 1]
      else 
        Qnew, hh = Two_Sum(Q, fnow)
        fnow = f[findex += 1]
	  end
      Q = Qnew
      if hh != 0.0
        h[hindex += 1] = hh
	  end
	end
  end

  while eindex <= elen
    Qnew, hh = Two_Sum(Q, enow)
    eindex += 1
    # We need an extra iteration to calculate Q
    # but we don't want to access e
    if eindex <= elen
        enow = e[eindex]
    end
    Q = Qnew
    if hh != 0.0
      h[hindex += 1] = hh
	end
  end

  while findex <= flen
    Qnew, hh = Two_Sum(Q, fnow)
    findex += 1
    if findex <= flen
        fnow = f[findex]
    end
    Q = Qnew
    if hh != 0.0
      h[hindex += 1] = hh
	end
  end
  if (Q != 0.0) || (hindex == 0)
    h[hindex += 1] = Q
  end
  return hindex
end
end