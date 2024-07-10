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
  @inline function fast_expansion_sum_zeroelim(elen::Int, e, flen::Int, f, ::Val{N})::Tuple{NTuple,Int} where {N}
    h = ntuple(i -> zero(typeof(e[1])), Val(N))
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
          #h[hindex+=1] = hh
          h = Base.setindex(h, hh, hindex+=1)
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
            #h[hindex += 1] = hh
            h = Base.setindex(h, hh, hindex+=1)
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
          #h[hindex += 1] = hh
          h = Base.setindex(h, hh, hindex+=1)
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
          #h[hindex += 1] = hh
          h = Base.setindex(h, hh, hindex+=1)
        end
      end
      if (Q != 0.0) || (hindex == 0)
        #h[hindex += 1] = Q
        h = Base.setindex(h, Q, hindex+=1)
      end
      return h, hindex
    end
  end