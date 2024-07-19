function inspherefast(pa, pb, pc, pd, pe)
    @inbounds begin
        aex = pa[1] - pe[1]
        bex = pb[1] - pe[1]
        cex = pc[1] - pe[1]
        dex = pd[1] - pe[1]
        aey = pa[2] - pe[2]
        bey = pb[2] - pe[2]
        cey = pc[2] - pe[2]
        dey = pd[2] - pe[2]
        aez = pa[3] - pe[3]
        bez = pb[3] - pe[3]
        cez = pc[3] - pe[3]
        dez = pd[3] - pe[3]

        ab = aex * bey - bex * aey
        bc = bex * cey - cex * bey
        cd = cex * dey - dex * cey
        da = dex * aey - aex * dey

        ac = aex * cey - cex * aey
        bd = bex * dey - dex * bey

        abc = aez * bc - bez * ac + cez * ab
        bcd = bez * cd - cez * bd + dez * bc
        cda = cez * da + dez * ac + aez * cd
        dab = dez * ab + aez * bd + bez * da

        alift = aex * aex + aey * aey + aez * aez
        blift = bex * bex + bey * bey + bez * bez
        clift = cex * cex + cey * cey + cez * cez
        dlift = dex * dex + dey * dey + dez * dez

        return (dlift * abc - clift * dab) + (blift * cda - alift * bcd)
    end
end

function insphereexact(pa, pb, pc, pd, pe)
    cache = InsphereCache{eltype(pa)}()
    return _insphereexact(pa, pb, pc, pd, pe, cache)
end
function _insphereexact(pa, pb, pc, pd, pe, cache)
    @inbounds begin
        axby1, axby0 = Two_Product(pa[1], pb[2])
        bxay1, bxay0 = Two_Product(pb[1], pa[2])
        ab3, ab2, ab1, ab0 = Two_Two_Diff(axby1, axby0, bxay1, bxay0)
        ab = (ab0, ab1, ab2, ab3)

        bxcy1, bxcy0 = Two_Product(pb[1], pc[2])
        cxby1, cxby0 = Two_Product(pc[1], pb[2])
        bc3, bc2, bc1, bc0 = Two_Two_Diff(bxcy1, bxcy0, cxby1, cxby0)
        bc = (bc0, bc1, bc2, bc3)

        cxdy1, cxdy0 = Two_Product(pc[1], pd[2])
        dxcy1, dxcy0 = Two_Product(pd[1], pc[2])
        cd3, cd2, cd1, cd0 = Two_Two_Diff(cxdy1, cxdy0, dxcy1, dxcy0)
        cd = (cd0, cd1, cd2, cd3)

        dxey1, dxey0 = Two_Product(pd[1], pe[2])
        exdy1, exdy0 = Two_Product(pe[1], pd[2])
        de3, de2, de1, de0 = Two_Two_Diff(dxey1, dxey0, exdy1, exdy0)
        de = (de0, de1, de2, de3)

        exay1, exay0 = Two_Product(pe[1], pa[2])
        axey1, axey0 = Two_Product(pa[1], pe[2])
        ea3, ea2, ea1, ea0 = Two_Two_Diff(exay1, exay0, axey1, axey0)
        ea = (ea0, ea1, ea2, ea3)

        axcy1, axcy0 = Two_Product(pa[1], pc[2])
        cxay1, cxay0 = Two_Product(pc[1], pa[2])
        ac3, ac2, ac1, ac0 = Two_Two_Diff(axcy1, axcy0, cxay1, cxay0)
        ac = (ac0, ac1, ac2, ac3)

        bxdy1, bxdy0 = Two_Product(pb[1], pd[2])
        dxby1, dxby0 = Two_Product(pd[1], pb[2])
        bd3, bd2, bd1, bd0 = Two_Two_Diff(bxdy1, bxdy0, dxby1, dxby0)
        bd = (bd0, bd1, bd2, bd3)

        cxey1, cxey0 = Two_Product(pc[1], pe[2])
        excy1, excy0 = Two_Product(pe[1], pc[2])
        ce3, ce2, ce1, ce0 = Two_Two_Diff(cxey1, cxey0, excy1, excy0)
        ce = (ce0, ce1, ce2, ce3)

        dxay1, dxay0 = Two_Product(pd[1], pa[2])
        axdy1, axdy0 = Two_Product(pa[1], pd[2])
        da3, da2, da1, da0 = Two_Two_Diff(dxay1, dxay0, axdy1, axdy0)
        da = (da0, da1, da2, da3)

        exby1, exby0 = Two_Product(pe[1], pb[2])
        bxey1, bxey0 = Two_Product(pb[1], pe[2])
        eb3, eb2, eb1, eb0 = Two_Two_Diff(exby1, exby0, bxey1, bxey0)
        eb = (eb0, eb1, eb2, eb3)

        temp8a, temp8alen = scale_expansion_zeroelim(4, bc, pa[3], cache.h8)
        temp8b, temp8blen = scale_expansion_zeroelim(4, ac, -pb[3], cache.h8)
        temp16, temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, cache.h16)
        temp8a, temp8alen = scale_expansion_zeroelim(4, ab, pc[3], temp8a)
        abc, abclen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, cache.h24)

        temp8a, temp8alen = scale_expansion_zeroelim(4, cd, pb[3], temp8a)
        temp8b, temp8blen = scale_expansion_zeroelim(4, bd, -pc[3], temp8b)
        temp16, temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16)
        temp8a, temp8alen = scale_expansion_zeroelim(4, bc, pd[3], temp8a)
        bcd, bcdlen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, cache.h24)

        temp8a, temp8alen = scale_expansion_zeroelim(4, de, pc[3], temp8a)
        temp8b, temp8blen = scale_expansion_zeroelim(4, ce, -pd[3], temp8b)
        temp16, temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16)
        temp8a, temp8alen = scale_expansion_zeroelim(4, cd, pe[3], temp8a)
        cde, cdelen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, cache.h24)

        temp8a, temp8alen = scale_expansion_zeroelim(4, ea, pd[3], temp8a)
        temp8b, temp8blen = scale_expansion_zeroelim(4, da, -pe[3], temp8b)
        temp16, temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16)
        temp8a, temp8alen = scale_expansion_zeroelim(4, de, pa[3], temp8a)
        dea, dealen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, cache.h24)

        temp8a, temp8alen = scale_expansion_zeroelim(4, ab, pe[3], temp8a)
        temp8b, temp8blen = scale_expansion_zeroelim(4, eb, -pa[3], temp8b)
        temp16, temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16)
        temp8a, temp8alen = scale_expansion_zeroelim(4, ea, pb[3], temp8a)
        eab, eablen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, cache.h24)

        temp8a, temp8alen = scale_expansion_zeroelim(4, bd, pa[3], temp8a)
        temp8b, temp8blen = scale_expansion_zeroelim(4, da, pb[3], temp8b)
        temp16, temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16)
        temp8a, temp8alen = scale_expansion_zeroelim(4, ab, pd[3], temp8a)
        abd, abdlen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, cache.h24)

        temp8a, temp8alen = scale_expansion_zeroelim(4, ce, pb[3], temp8a)
        temp8b, temp8blen = scale_expansion_zeroelim(4, eb, pc[3], temp8b)
        temp16, temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16)
        temp8a, temp8alen = scale_expansion_zeroelim(4, bc, pe[3], temp8a)
        bce, bcelen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, cache.h24)

        temp8a, temp8alen = scale_expansion_zeroelim(4, da, pc[3], temp8a)
        temp8b, temp8blen = scale_expansion_zeroelim(4, ac, pd[3], temp8b)
        temp16, temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16)
        temp8a, temp8alen = scale_expansion_zeroelim(4, cd, pa[3], temp8a)
        cda, cdalen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, cache.h24)

        temp8a, temp8alen = scale_expansion_zeroelim(4, eb, pd[3], temp8a)
        temp8b, temp8blen = scale_expansion_zeroelim(4, bd, pe[3], temp8b)
        temp16, temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16)
        temp8a, temp8alen = scale_expansion_zeroelim(4, de, pb[3], temp8a)
        deb, deblen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, cache.h24)

        temp8a, temp8alen = scale_expansion_zeroelim(4, ac, pe[3], temp8a)
        temp8b, temp8blen = scale_expansion_zeroelim(4, ce, pa[3], temp8b)
        temp16, temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16)
        temp8a, temp8alen = scale_expansion_zeroelim(4, ea, pc[3], temp8a)
        eac, eaclen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, cache.h24)

        temp48a, temp48alen = fast_expansion_sum_zeroelim(cdelen, cde, bcelen, bce, cache.h48_1)
        temp48b, temp48blen = fast_expansion_sum_zeroelim(deblen, deb, bcdlen, bcd, cache.h48_2)
        for i in 1:temp48blen
            temp48b[i] = -temp48b[i]
        end
        bcde, bcdelen = fast_expansion_sum_zeroelim(temp48alen, temp48a, temp48blen, temp48b, cache.h96_1)
        temp192, xlen = scale_expansion_zeroelim(bcdelen, bcde, pa[1], cache.h192)
        det384x, xlen = scale_expansion_zeroelim(xlen, temp192, pa[1], cache.h384_1)
        temp192, ylen = scale_expansion_zeroelim(bcdelen, bcde, pa[2], temp192)
        det384y, ylen = scale_expansion_zeroelim(ylen, temp192, pa[2], cache.h384_2)
        temp192, zlen = scale_expansion_zeroelim(bcdelen, bcde, pa[3], temp192)
        det384z, zlen = scale_expansion_zeroelim(zlen, temp192, pa[3], cache.h384_3)
        detxy, xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, cache.h768_1)
        adet, alen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, cache.h1152_1)

        temp48a, temp48alen = fast_expansion_sum_zeroelim(dealen, dea, cdalen, cda, temp48a)
        temp48b, temp48blen = fast_expansion_sum_zeroelim(eaclen, eac, cdelen, cde, temp48b)
        for i in 1:temp48blen
            temp48b[i] = -temp48b[i]
        end
        cdea, cdealen = fast_expansion_sum_zeroelim(temp48alen, temp48a, temp48blen, temp48b, cache.h96_2)
        temp192, xlen = scale_expansion_zeroelim(cdealen, cdea, pb[1], temp192)
        det384x, xlen = scale_expansion_zeroelim(xlen, temp192, pb[1], det384x)
        temp192, ylen = scale_expansion_zeroelim(cdealen, cdea, pb[2], temp192)
        det384y, ylen = scale_expansion_zeroelim(ylen, temp192, pb[2], det384y)
        temp192, zlen = scale_expansion_zeroelim(cdealen, cdea, pb[3], temp192)
        det384z, zlen = scale_expansion_zeroelim(zlen, temp192, pb[3], det384z)
        detxy, xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy)
        bdet, blen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, cache.h1152_2)

        temp48a, temp48alen = fast_expansion_sum_zeroelim(eablen, eab, deblen, deb, temp48a)
        temp48b, temp48blen = fast_expansion_sum_zeroelim(abdlen, abd, dealen, dea, temp48b)
        for i in 1:temp48blen
            temp48b[i] = -temp48b[i]
        end
        deab, deablen = fast_expansion_sum_zeroelim(temp48alen, temp48a, temp48blen, temp48b, cache.h96_3)
        temp192, xlen = scale_expansion_zeroelim(deablen, deab, pc[1], temp192)
        det384x, xlen = scale_expansion_zeroelim(xlen, temp192, pc[1], det384x)
        temp192, ylen = scale_expansion_zeroelim(deablen, deab, pc[2], temp192)
        det384y, ylen = scale_expansion_zeroelim(ylen, temp192, pc[2], det384y)
        temp192, zlen = scale_expansion_zeroelim(deablen, deab, pc[3], temp192)
        det384z, zlen = scale_expansion_zeroelim(zlen, temp192, pc[3], det384z)
        detxy, xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy)
        cdet, clen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, cache.h1152_3)

        temp48a, temp48alen = fast_expansion_sum_zeroelim(abclen, abc, eaclen, eac, temp48a)
        temp48b, temp48blen = fast_expansion_sum_zeroelim(bcelen, bce, eablen, eab, temp48b)
        for i in 1:temp48blen
            temp48b[i] = -temp48b[i]
        end
        eabc, eabclen = fast_expansion_sum_zeroelim(temp48alen, temp48a, temp48blen, temp48b, cache.h96_4)
        temp192, xlen = scale_expansion_zeroelim(eabclen, eabc, pd[1], temp192)
        det384x, xlen = scale_expansion_zeroelim(xlen, temp192, pd[1], det384x)
        temp192, ylen = scale_expansion_zeroelim(eabclen, eabc, pd[2], temp192)
        det384y, ylen, = scale_expansion_zeroelim(ylen, temp192, pd[2], det384y)
        temp192, zlen = scale_expansion_zeroelim(eabclen, eabc, pd[3], temp192)
        det384z, zlen = scale_expansion_zeroelim(zlen, temp192, pd[3], det384z)
        detxy, xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy)
        ddet, dlen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, cache.h1152_4)

        temp48a, temp48alen = fast_expansion_sum_zeroelim(bcdlen, bcd, abdlen, abd, temp48a)
        temp48b, temp48blen = fast_expansion_sum_zeroelim(cdalen, cda, abclen, abc, temp48b)
        for i in 1:temp48blen
            temp48b[i] = -temp48b[i]
        end
        abcd, abcdlen = fast_expansion_sum_zeroelim(temp48alen, temp48a, temp48blen, temp48b, cache.h96_5)
        temp192, xlen = scale_expansion_zeroelim(abcdlen, abcd, pe[1], temp192)
        det384x, xlen = scale_expansion_zeroelim(xlen, temp192, pe[1], det384x)
        temp192, ylen = scale_expansion_zeroelim(abcdlen, abcd, pe[2], temp192)
        det384y, ylen = scale_expansion_zeroelim(ylen, temp192, pe[2], det384y)
        temp192, zlen = scale_expansion_zeroelim(abcdlen, abcd, pe[3], temp192)
        det384z, zlen = scale_expansion_zeroelim(zlen, temp192, pe[3], det384z)
        detxy, xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy)
        edet, elen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, cache.h1152_5)

        abdet, ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, cache.h2304_1)
        cddet, cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cache.h2304_2)
        cdedet, cdelen = fast_expansion_sum_zeroelim(cdlen, cddet, elen, edet, cache.h3456)
        deter, deterlen = fast_expansion_sum_zeroelim(ablen, abdet, cdelen, cdedet, cache.h5760)

        return deter[deterlen]
    end
end

function insphereslow(pa, pb, pc, pd, pe)
    cache = InsphereCache{eltype(pa)}()
    return _insphereslow(pa, pb, pc, pd, pe, cache)
end
function _insphereslow(pa, pb, pc, pd, pe, cache)
    @inbounds begin
        aex, aextail = Two_Diff(pa[1], pe[1])
        aey, aeytail = Two_Diff(pa[2], pe[2])
        aez, aeztail = Two_Diff(pa[3], pe[3])
        bex, bextail = Two_Diff(pb[1], pe[1])
        bey, beytail = Two_Diff(pb[2], pe[2])
        bez, beztail = Two_Diff(pb[3], pe[3])
        cex, cextail = Two_Diff(pc[1], pe[1])
        cey, ceytail = Two_Diff(pc[2], pe[2])
        cez, ceztail = Two_Diff(pc[3], pe[3])
        dex, dextail = Two_Diff(pd[1], pe[1])
        dey, deytail = Two_Diff(pd[2], pe[2])
        dez, deztail = Two_Diff(pd[3], pe[3])

        axby7, axby6, axby5, axby4, axby3, axby2, axby1, axby0 = Two_Two_Product(aex, aextail, bey, beytail)
        axby = (axby0, axby1, axby2, axby3, axby4, axby5, axby6, axby7)
        negate = -aey
        negatetail = -aeytail
        bxay7, bxay6, bxay5, bxay4, bxay3, bxay2, bxay1, bxay0 = Two_Two_Product(bex, bextail, negate, negatetail)
        bxay = (bxay0, bxay1, bxay2, bxay3, bxay4, bxay5, bxay6, bxay7)
        ab, ablen = fast_expansion_sum_zeroelim(8, axby, 8, bxay, cache.h16)
        bxcy7, bxcy6, bxcy5, bxcy4, bxcy3, bxcy2, bxcy1, bxcy0 = Two_Two_Product(bex, bextail, cey, ceytail)
        bxcy = (bxcy0, bxcy1, bxcy2, bxcy3, bxcy4, bxcy5, bxcy6, bxcy7)
        negate = -bey
        negatetail = -beytail
        cxby7, cxby6, cxby5, cxby4, cxby3, cxby2, cxby1, cxby0 = Two_Two_Product(cex, cextail, negate, negatetail)
        cxby = (cxby0, cxby1, cxby2, cxby3, cxby4, cxby5, cxby6, cxby7)
        bc, bclen = fast_expansion_sum_zeroelim(8, bxcy, 8, cxby, cache.h16)
        cxdy7, cxdy6, cxdy5, cxdy4, cxdy3, cxdy2, cxdy1, cxdy0 = Two_Two_Product(cex, cextail, dey, deytail)
        cxdy = (cxdy0, cxdy1, cxdy2, cxdy3, cxdy4, cxdy5, cxdy6, cxdy7)
        negate = -cey
        negatetail = -ceytail
        dxcy7, dxcy6, dxcy5, dxcy4, dxcy3, dxcy2, dxcy1, dxcy0 = Two_Two_Product(dex, dextail, negate, negatetail)
        dxcy = (dxcy0, dxcy1, dxcy2, dxcy3, dxcy4, dxcy5, dxcy6, dxcy7)
        cd, cdlen = fast_expansion_sum_zeroelim(8, cxdy, 8, dxcy, cache.h16)
        dxay7, dxay6, dxay5, dxay4, dxay3, dxay2, dxay1, dxay0 = Two_Two_Product(dex, dextail, aey, aeytail)
        dxay = (dxay0, dxay1, dxay2, dxay3, dxay4, dxay5, dxay6, dxay7)
        negate = -dey
        negatetail = -deytail
        axdy7, axdy6, axdy5, axdy4, axdy3, axdy2, axdy1, axdy0 = Two_Two_Product(aex, aextail, negate, negatetail)
        axdy = (axdy0, axdy1, axdy2, axdy3, axdy4, axdy5, axdy6, axdy7)
        da, dalen = fast_expansion_sum_zeroelim(8, dxay, 8, axdy, cache.h16)
        axcy7, axcy6, axcy5, axcy4, axcy3, axcy2, axcy1, axcy0 = Two_Two_Product(aex, aextail, cey, ceytail)
        axcy = (axcy0, axcy1, axcy2, axcy3, axcy4, axcy5, axcy6, axcy7)
        negate = -aey
        negatetail = -aeytail
        cxay7, cxay6, cxay5, cxay4, cxay3, cxay2, cxay1, cxay0 = Two_Two_Product(cex, cextail, negate, negatetail)
        cxay = (cxay0, cxay1, cxay2, cxay3, cxay4, cxay5, cxay6, cxay7)
        ac, aclen = fast_expansion_sum_zeroelim(8, axcy, 8, cxay, cache.h16)
        bxdy7, bxdy6, bxdy5, bxdy4, bxdy3, bxdy2, bxdy1, bxdy0 = Two_Two_Product(bex, bextail, dey, deytail)
        bxdy = (bxdy0, bxdy1, bxdy2, bxdy3, bxdy4, bxdy5, bxdy6, bxdy7)
        negate = -bey
        negatetail = -beytail
        dxby7, dxby6, dxby5, dxby4, dxby3, dxby2, dxby1, dxby0 = Two_Two_Product(dex, dextail, negate, negatetail)
        dxby = (dxby0, dxby1, dxby2, dxby3, dxby4, dxby5, dxby6, dxby7)
        bd, bdlen = fast_expansion_sum_zeroelim(8, bxdy, 8, dxby, cache.h16)

        temp32a, temp32alen = scale_expansion_zeroelim(cdlen, cd, -bez, cache.h32)
        temp32b, temp32blen = scale_expansion_zeroelim(cdlen, cd, -beztail, cache.h32)
        temp64a, temp64alen = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, cache.h64_1)
        temp32a, temp32alen = scale_expansion_zeroelim(bdlen, bd, cez, temp32a)
        temp32b, temp32blen = scale_expansion_zeroelim(bdlen, bd, ceztail, temp32b)
        temp64b, temp64blen = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, cache.h64_2)
        temp32a, temp32alen = scale_expansion_zeroelim(bclen, bc, -dez, temp32a)
        temp32b, temp32blen = scale_expansion_zeroelim(bclen, bc, -deztail, temp32b)
        temp64c, temp64clen = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, cache.h64_3)
        temp128, temp128len = fast_expansion_sum_zeroelim(temp64alen, temp64a, temp64blen, temp64b, cache.h128)
        temp192, temp192len = fast_expansion_sum_zeroelim(temp64clen, temp64c, temp128len, temp128, cache.h192)
        detx, xlen = scale_expansion_zeroelim(temp192len, temp192, aex, cache.h384_1)
        detxx, xxlen = scale_expansion_zeroelim(xlen, detx, aex, cache.h768_1)
        detxt, xtlen = scale_expansion_zeroelim(temp192len, temp192, aextail, cache.h384_2)
        detxxt, xxtlen = scale_expansion_zeroelim(xtlen, detxt, aex, cache.h768_2)
        for i in 1:xxtlen
            detxxt[i] *= 2.0
        end
        detxtxt, xtxtlen = scale_expansion_zeroelim(xtlen, detxt, aextail, cache.h768_3)
        x1, x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, cache.h1536_1)
        x2, x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, cache.h2304_1)
        dety, ylen = scale_expansion_zeroelim(temp192len, temp192, aey, cache.h384_3)
        detyy, yylen = scale_expansion_zeroelim(ylen, dety, aey, cache.h768_4)
        detyt, ytlen = scale_expansion_zeroelim(temp192len, temp192, aeytail, cache.h384_4)
        detyyt, yytlen = scale_expansion_zeroelim(ytlen, detyt, aey, cache.h768_5)
        for i in 1:yytlen
            detyyt[i] *= 2.0
        end
        detytyt, ytytlen = scale_expansion_zeroelim(ytlen, detyt, aeytail, cache.h768_6)
        y1, y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, cache.h1536_2)
        y2, y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, cache.h2304_2)
        detz, zlen = scale_expansion_zeroelim(temp192len, temp192, aez, cache.h384_5)
        detzz, zzlen = scale_expansion_zeroelim(zlen, detz, aez, cache.h768_7)
        detzt, ztlen = scale_expansion_zeroelim(temp192len, temp192, aeztail, cache.h384_6)
        detzzt, zztlen = scale_expansion_zeroelim(ztlen, detzt, aez, cache.h768_8)
        for i in 1:zztlen
            detzzt[i] *= 2.0
        end
        detztzt, ztztlen = scale_expansion_zeroelim(ztlen, detzt, aeztail, cache.h768_9)
        z1, z1len = fast_expansion_sum_zeroelim(zzlen, detzz, zztlen, detzzt, cache.h1536_3)
        z2, z2len = fast_expansion_sum_zeroelim(z1len, z1, ztztlen, detztzt, cache.h2304_3)
        detxy, xylen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, cache.h4608)
        adet, alen = fast_expansion_sum_zeroelim(z2len, z2, xylen, detxy, cache.h6912_1)

        temp32a, temp32alen = scale_expansion_zeroelim(dalen, da, cez, temp32a)
        temp32b, temp32blen = scale_expansion_zeroelim(dalen, da, ceztail, temp32b)
        temp64a, temp64alen = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, temp64a)
        temp32a, temp32alen = scale_expansion_zeroelim(aclen, ac, dez, temp32a)
        temp32b, temp32blen = scale_expansion_zeroelim(aclen, ac, deztail, temp32b)
        temp64b, temp64blen = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, temp64b)
        temp32a, temp32alen = scale_expansion_zeroelim(cdlen, cd, aez, temp32a)
        temp32b, temp32blen = scale_expansion_zeroelim(cdlen, cd, aeztail, temp32b)
        temp64c, temp64clen = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, temp64c)
        temp128, temp128len = fast_expansion_sum_zeroelim(temp64alen, temp64a, temp64blen, temp64b, temp128)
        temp192, temp192len = fast_expansion_sum_zeroelim(temp64clen, temp64c, temp128len, temp128, temp192)
        detx, xlen = scale_expansion_zeroelim(temp192len, temp192, bex, detx)
        detxx, xxlen = scale_expansion_zeroelim(xlen, detx, bex, detxx)
        detxt, xtlen = scale_expansion_zeroelim(temp192len, temp192, bextail, detxt)
        detxxt, xxtlen = scale_expansion_zeroelim(xtlen, detxt, bex, detxxt)
        for i in 1:xxtlen
            detxxt[i] *= 2.0
        end
        detxtxt, xtxtlen = scale_expansion_zeroelim(xtlen, detxt, bextail, detxtxt)
        x1, x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, x1)
        x2, x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, x2)
        dety, ylen = scale_expansion_zeroelim(temp192len, temp192, bey, dety)
        detyy, yylen = scale_expansion_zeroelim(ylen, dety, bey, detyy)
        detyt, ytlen = scale_expansion_zeroelim(temp192len, temp192, beytail, detyt)
        detyyt, yytlen = scale_expansion_zeroelim(ytlen, detyt, bey, detyyt)
        for i in 1:yytlen
            detyyt[i] *= 2.0
        end
        detytyt, ytytlen = scale_expansion_zeroelim(ytlen, detyt, beytail, detytyt)
        y1, y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, y1)
        y2, y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, y2)
        detz, zlen = scale_expansion_zeroelim(temp192len, temp192, bez, detz)
        detzz, zzlen = scale_expansion_zeroelim(zlen, detz, bez, detzz)
        detzt, ztlen = scale_expansion_zeroelim(temp192len, temp192, beztail, detzt)
        detzzt, zztlen = scale_expansion_zeroelim(ztlen, detzt, bez, detzzt)
        for i in 1:zztlen
            detzzt[i] *= 2.0
        end
        detztzt, ztztlen = scale_expansion_zeroelim(ztlen, detzt, beztail, detztzt)
        z1, z1len = fast_expansion_sum_zeroelim(zzlen, detzz, zztlen, detzzt, z1)
        z2, z2len = fast_expansion_sum_zeroelim(z1len, z1, ztztlen, detztzt, z2)
        detxy, xylen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, detxy)
        bdet, blen = fast_expansion_sum_zeroelim(z2len, z2, xylen, detxy, cache.h6912_2)

        temp32a, temp32alen = scale_expansion_zeroelim(ablen, ab, -dez, temp32a)
        temp32b, temp32blen = scale_expansion_zeroelim(ablen, ab, -deztail, temp32b)
        temp64a, temp64alen = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, temp64a)
        temp32a, temp32alen = scale_expansion_zeroelim(bdlen, bd, -aez, temp32a)
        temp32b, temp32blen = scale_expansion_zeroelim(bdlen, bd, -aeztail, temp32b)
        temp64b, temp64blen = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, temp64b)
        temp32a, temp32alen = scale_expansion_zeroelim(dalen, da, -bez, temp32a)
        temp32b, temp32blen = scale_expansion_zeroelim(dalen, da, -beztail, temp32b)
        temp64c, temp64clen = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, temp64c)
        temp128, temp128len = fast_expansion_sum_zeroelim(temp64alen, temp64a, temp64blen, temp64b, temp128)
        temp192, temp192len = fast_expansion_sum_zeroelim(temp64clen, temp64c, temp128len, temp128, temp192)
        detx, xlen = scale_expansion_zeroelim(temp192len, temp192, cex, detx)
        detxx, xxlen = scale_expansion_zeroelim(xlen, detx, cex, detxx)
        detxt, xtlen = scale_expansion_zeroelim(temp192len, temp192, cextail, detxt)
        detxxt, xxtlen = scale_expansion_zeroelim(xtlen, detxt, cex, detxxt)
        for i in 1:xxtlen
            detxxt[i] *= 2.0
        end
        detxtxt, xtxtlen = scale_expansion_zeroelim(xtlen, detxt, cextail, detxtxt)
        x1, x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, x1)
        x2, x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, x2)
        dety, ylen = scale_expansion_zeroelim(temp192len, temp192, cey, dety)
        detyy, yylen = scale_expansion_zeroelim(ylen, dety, cey, detyy)
        detyt, ytlen = scale_expansion_zeroelim(temp192len, temp192, ceytail, detyt)
        detyyt, yytlen = scale_expansion_zeroelim(ytlen, detyt, cey, detyyt)
        for i in 1:yytlen
            detyyt[i] *= 2.0
        end
        detytyt, ytytlen = scale_expansion_zeroelim(ytlen, detyt, ceytail, detytyt)
        y1, y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, y1)
        y2, y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, y2)
        detz, zlen = scale_expansion_zeroelim(temp192len, temp192, cez, detz)
        detzz, zzlen = scale_expansion_zeroelim(zlen, detz, cez, detzz)
        detzt, ztlen = scale_expansion_zeroelim(temp192len, temp192, ceztail, detzt)
        detzzt, zztlen = scale_expansion_zeroelim(ztlen, detzt, cez, detzzt)
        for i in 1:zztlen
            detzzt[i] *= 2.0
        end
        detztzt, ztztlen = scale_expansion_zeroelim(ztlen, detzt, ceztail, detztzt)
        z1, z1len = fast_expansion_sum_zeroelim(zzlen, detzz, zztlen, detzzt, z1)
        z2, z2len = fast_expansion_sum_zeroelim(z1len, z1, ztztlen, detztzt, z2)
        detxy, xylen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, detxy)
        cdet, clen = fast_expansion_sum_zeroelim(z2len, z2, xylen, detxy, cache.h6912_3)

        temp32a, temp32alen = scale_expansion_zeroelim(bclen, bc, aez, temp32a)
        temp32b, temp32blen = scale_expansion_zeroelim(bclen, bc, aeztail, temp32b)
        temp64a, temp64alen = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, temp64a)
        temp32a, temp32alen = scale_expansion_zeroelim(aclen, ac, -bez, temp32a)
        temp32b, temp32blen = scale_expansion_zeroelim(aclen, ac, -beztail, temp32b)
        temp64b, temp64blen = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, temp64b)
        temp32a, temp32alen = scale_expansion_zeroelim(ablen, ab, cez, temp32a)
        temp32b, temp32blen = scale_expansion_zeroelim(ablen, ab, ceztail, temp32b)
        temp64c, temp64clen = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, temp64c)
        temp128, temp128len = fast_expansion_sum_zeroelim(temp64alen, temp64a, temp64blen, temp64b, temp128)
        temp192, temp192len = fast_expansion_sum_zeroelim(temp64clen, temp64c, temp128len, temp128, temp192)
        detx, xlen = scale_expansion_zeroelim(temp192len, temp192, dex, detx)
        detxx, xxlen = scale_expansion_zeroelim(xlen, detx, dex, detxx)
        detxt, xtlen = scale_expansion_zeroelim(temp192len, temp192, dextail, detxt)
        detxxt, xxtlen = scale_expansion_zeroelim(xtlen, detxt, dex, detxxt)
        for i in 1:xxtlen
            detxxt[i] *= 2.0
        end
        detxtxt, xtxtlen = scale_expansion_zeroelim(xtlen, detxt, dextail, detxtxt)
        x1, x1len = fast_expansion_sum_zeroelim(xxlen, detxx, xxtlen, detxxt, x1)
        x2, x2len = fast_expansion_sum_zeroelim(x1len, x1, xtxtlen, detxtxt, x2)
        dety, ylen = scale_expansion_zeroelim(temp192len, temp192, dey, dety)
        detyy, yylen = scale_expansion_zeroelim(ylen, dety, dey, detyy)
        detyt, ytlen = scale_expansion_zeroelim(temp192len, temp192, deytail, detyt)
        detyyt, yytlen = scale_expansion_zeroelim(ytlen, detyt, dey, detyyt)
        for i in 1:yytlen
            detyyt[i] *= 2.0
        end
        detytyt, ytytlen = scale_expansion_zeroelim(ytlen, detyt, deytail, detytyt)
        y1, y1len = fast_expansion_sum_zeroelim(yylen, detyy, yytlen, detyyt, y1)
        y2, y2len = fast_expansion_sum_zeroelim(y1len, y1, ytytlen, detytyt, y2)
        detz, zlen = scale_expansion_zeroelim(temp192len, temp192, dez, detz)
        detzz, zzlen = scale_expansion_zeroelim(zlen, detz, dez, detzz)
        detzt, ztlen = scale_expansion_zeroelim(temp192len, temp192, deztail, detzt)
        detzzt, zztlen = scale_expansion_zeroelim(ztlen, detzt, dez, detzzt)
        for i in 1:zztlen
            detzzt[i] *= 2.0
        end
        detztzt, ztztlen = scale_expansion_zeroelim(ztlen, detzt, deztail, detztzt)
        z1, z1len = fast_expansion_sum_zeroelim(zzlen, detzz, zztlen, detzzt, z1)
        z2, z2len = fast_expansion_sum_zeroelim(z1len, z1, ztztlen, detztzt, z2)
        detxy, xylen = fast_expansion_sum_zeroelim(x2len, x2, y2len, y2, detxy)
        ddet, dlen = fast_expansion_sum_zeroelim(z2len, z2, xylen, detxy, cache.h6912_4)

        abdet, ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, cache.h13824_1)
        cddet, cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cache.h13824_2)
        deter, deterlen = fast_expansion_sum_zeroelim(ablen, abdet, cdlen, cddet, cache.h27648)

        return deter[deterlen]
    end
end

function insphereadapt(pa, pb, pc, pd, pe, permanent)
    cache = InsphereCache{eltype(pa)}()
    return _insphereadapt(pa, pb, pc, pd, pe, permanent, cache)
end
function _insphereadapt(pa, pb, pc, pd, pe, permanent, cache::InsphereCache{T}) where {T}
    @inbounds begin
        aex = pa[1] - pe[1]
        bex = pb[1] - pe[1]
        cex = pc[1] - pe[1]
        dex = pd[1] - pe[1]
        aey = pa[2] - pe[2]
        bey = pb[2] - pe[2]
        cey = pc[2] - pe[2]
        dey = pd[2] - pe[2]
        aez = pa[3] - pe[3]
        bez = pb[3] - pe[3]
        cez = pc[3] - pe[3]
        dez = pd[3] - pe[3]

        aexbey1, aexbey0 = Two_Product(aex, bey)
        bexaey1, bexaey0 = Two_Product(bex, aey)
        ab3, ab2, ab1, ab0 = Two_Two_Diff(aexbey1, aexbey0, bexaey1, bexaey0)
        ab = (ab0, ab1, ab2, ab3)

        bexcey1, bexcey0 = Two_Product(bex, cey)
        cexbey1, cexbey0 = Two_Product(cex, bey)
        bc3, bc2, bc1, bc0 = Two_Two_Diff(bexcey1, bexcey0, cexbey1, cexbey0)
        bc = (bc0, bc1, bc2, bc3)

        cexdey1, cexdey0 = Two_Product(cex, dey)
        dexcey1, dexcey0 = Two_Product(dex, cey)
        cd3, cd2, cd1, cd0 = Two_Two_Diff(cexdey1, cexdey0, dexcey1, dexcey0)
        cd = (cd0, cd1, cd2, cd3)

        dexaey1, dexaey0 = Two_Product(dex, aey)
        aexdey1, aexdey0 = Two_Product(aex, dey)
        da3, da2, da1, da0 = Two_Two_Diff(dexaey1, dexaey0, aexdey1, aexdey0)
        da = (da0, da1, da2, da3)

        aexcey1, aexcey0 = Two_Product(aex, cey)
        cexaey1, cexaey0 = Two_Product(cex, aey)
        ac3, ac2, ac1, ac0 = Two_Two_Diff(aexcey1, aexcey0, cexaey1, cexaey0)
        ac = (ac0, ac1, ac2, ac3)

        bexdey1, bexdey0 = Two_Product(bex, dey)
        dexbey1, dexbey0 = Two_Product(dex, bey)
        bd3, bd2, bd1, bd0 = Two_Two_Diff(bexdey1, bexdey0, dexbey1, dexbey0)
        bd = (bd0, bd1, bd2, bd3)

        temp8a, temp8alen = scale_expansion_zeroelim(4, cd, bez, cache.h8)
        temp8b, temp8blen = scale_expansion_zeroelim(4, bd, -cez, cache.h8)
        temp8c, temp8clen = scale_expansion_zeroelim(4, bc, dez, cache.h8)
        temp16, temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, cache.h16)
        temp24, temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c, temp16len, temp16, cache.h24)
        temp48, temp48len = scale_expansion_zeroelim(temp24len, temp24, aex, cache.h48_1)
        xdet, xlen = scale_expansion_zeroelim(temp48len, temp48, -aex, cache.h96_1)
        temp48, temp48len = scale_expansion_zeroelim(temp24len, temp24, aey, temp48)
        ydet, ylen = scale_expansion_zeroelim(temp48len, temp48, -aey, cache.h96_2)
        temp48, temp48len = scale_expansion_zeroelim(temp24len, temp24, aez, temp48)
        zdet, zlen = scale_expansion_zeroelim(temp48len, temp48, -aez, cache.h96_3)
        xydet, xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, cache.h192)
        adet, alen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, cache.h288_1)

        temp8a, temp8alen = scale_expansion_zeroelim(4, da, cez, temp8a)
        temp8b, temp8blen = scale_expansion_zeroelim(4, ac, dez, temp8b)
        temp8c, temp8clen = scale_expansion_zeroelim(4, cd, aez, temp8c)
        temp16, temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16)
        temp24, temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c, temp16len, temp16, temp24)
        temp48, temp48len = scale_expansion_zeroelim(temp24len, temp24, bex, temp48)
        xdet, xlen = scale_expansion_zeroelim(temp48len, temp48, bex, xdet)
        temp48, temp48len = scale_expansion_zeroelim(temp24len, temp24, bey, temp48)
        ydet, ylen = scale_expansion_zeroelim(temp48len, temp48, bey, ydet)
        temp48, temp48len = scale_expansion_zeroelim(temp24len, temp24, bez, temp48)
        zdet, zlen = scale_expansion_zeroelim(temp48len, temp48, bez, zdet)
        xydet, xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, xydet)
        bdet, blen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, cache.h288_2)

        temp8a, temp8alen = scale_expansion_zeroelim(4, ab, dez, temp8a)
        temp8b, temp8blen = scale_expansion_zeroelim(4, bd, aez, temp8b)
        temp8c, temp8clen = scale_expansion_zeroelim(4, da, bez, temp8c)
        temp16, temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16)
        temp24, temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c, temp16len, temp16, temp24)
        temp48, temp48len = scale_expansion_zeroelim(temp24len, temp24, cex, temp48)
        xdet, xlen = scale_expansion_zeroelim(temp48len, temp48, -cex, xdet)
        temp48, temp48len = scale_expansion_zeroelim(temp24len, temp24, cey, temp48)
        ydet, ylen = scale_expansion_zeroelim(temp48len, temp48, -cey, ydet)
        temp48, temp48len = scale_expansion_zeroelim(temp24len, temp24, cez, temp48)
        zdet, zlen = scale_expansion_zeroelim(temp48len, temp48, -cez, zdet)
        xydet, xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, xydet)
        cdet, clen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, cache.h288_3)

        temp8a, temp8alen = scale_expansion_zeroelim(4, bc, aez, temp8a)
        temp8b, temp8blen = scale_expansion_zeroelim(4, ac, -bez, temp8b)
        temp8c, temp8clen = scale_expansion_zeroelim(4, ab, cez, temp8c)
        temp16, temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16)
        temp24, temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c, temp16len, temp16, temp24)
        temp48, temp48len = scale_expansion_zeroelim(temp24len, temp24, dex, temp48)
        xdet, xlen = scale_expansion_zeroelim(temp48len, temp48, dex, xdet)
        temp48, temp48len = scale_expansion_zeroelim(temp24len, temp24, dey, temp48)
        ydet, ylen = scale_expansion_zeroelim(temp48len, temp48, dey, ydet)
        temp48, temp48len = scale_expansion_zeroelim(temp24len, temp24, dez, temp48)
        zdet, zlen = scale_expansion_zeroelim(temp48len, temp48, dez, zdet)
        xydet, xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, xydet)
        ddet, dlen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, cache.h288_4)

        abdet, ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, cache.h576_1)
        cddet, cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cache.h576_2)
        fin1, finlength = fast_expansion_sum_zeroelim(ablen, abdet, cdlen, cddet, cache.h1152_1)

        det = estimate(finlength, fin1)
        errbound = isperrboundB(T) * permanent
        if (det ≥ errbound) || (-det ≥ errbound)
            return det
        end

        aextail = Two_Diff_Tail(pa[1], pe[1], aex)
        aeytail = Two_Diff_Tail(pa[2], pe[2], aey)
        aeztail = Two_Diff_Tail(pa[3], pe[3], aez)
        bextail = Two_Diff_Tail(pb[1], pe[1], bex)
        beytail = Two_Diff_Tail(pb[2], pe[2], bey)
        beztail = Two_Diff_Tail(pb[3], pe[3], bez)
        cextail = Two_Diff_Tail(pc[1], pe[1], cex)
        ceytail = Two_Diff_Tail(pc[2], pe[2], cey)
        ceztail = Two_Diff_Tail(pc[3], pe[3], cez)
        dextail = Two_Diff_Tail(pd[1], pe[1], dex)
        deytail = Two_Diff_Tail(pd[2], pe[2], dey)
        deztail = Two_Diff_Tail(pd[3], pe[3], dez)
        if iszero(aextail) && iszero(aeytail) && iszero(aeztail) &&
           iszero(bextail) && iszero(beytail) && iszero(beztail) &&
           iszero(cextail) && iszero(ceytail) && iszero(ceztail) &&
           iszero(dextail) && iszero(deytail) && iszero(deztail)
            return det
        end

        errbound = isperrboundC(T) * permanent + resulterrbound(T) * Absolute(det)
        abeps = (aex * beytail + bey * aextail) -
                (aey * bextail + bex * aeytail)
        bceps = (bex * ceytail + cey * bextail) -
                (bey * cextail + cex * beytail)
        cdeps = (cex * deytail + dey * cextail) -
                (cey * dextail + dex * ceytail)
        daeps = (dex * aeytail + aey * dextail) -
                (dey * aextail + aex * deytail)
        aceps = (aex * ceytail + cey * aextail) -
                (aey * cextail + cex * aeytail)
        bdeps = (bex * deytail + dey * bextail) -
                (bey * dextail + dex * beytail)
        detadd = (((bex * bex + bey * bey + bez * bez) *
                   ((cez * daeps + dez * aceps + aez * cdeps) +
                    (ceztail * da3 + deztail * ac3 + aeztail * cd3)) +
                   (dex * dex + dey * dey + dez * dez) *
                   ((aez * bceps - bez * aceps + cez * abeps) +
                    (aeztail * bc3 - beztail * ac3 + ceztail * ab3))) -
                  ((aex * aex + aey * aey + aez * aez) *
                   ((bez * cdeps - cez * bdeps + dez * bceps) +
                    (beztail * cd3 - ceztail * bd3 + deztail * bc3)) +
                   (cex * cex + cey * cey + cez * cez) *
                   ((dez * abeps + aez * bdeps + bez * daeps) +
                    (deztail * ab3 + aeztail * bd3 + beztail * da3)))) +
                 2.0 * (((bex * bextail + bey * beytail + bez * beztail) *
                         (cez * da3 + dez * ac3 + aez * cd3) +
                         (dex * dextail + dey * deytail + dez * deztail) *
                         (aez * bc3 - bez * ac3 + cez * ab3)) -
                        ((aex * aextail + aey * aeytail + aez * aeztail) *
                         (bez * cd3 - cez * bd3 + dez * bc3) +
                         (cex * cextail + cey * ceytail + cez * ceztail) *
                         (dez * ab3 + aez * bd3 + bez * da3)))
        det = T(det + detadd) # Had to change this to match how C handles the 2.0 multiplication with Float32
        if (det ≥ errbound) || (-det ≥ errbound)
            return det
        end

        return _insphereexact(pa, pb, pc, pd, pe, cache)
        #=
        Interesting note: Shewchuk's paper on p. 352 says
            This implementation differs from the other
            tests in that, due to programmer laziness, D is not computed incrementally from B;
            rather, if C is not accurate enough, D is computed from scratch. Fortunately, C is usually
            accurate enough.    
        This is why this is the only one of the four predicates that just uses the exact computation 
        when the result is not known yet - laziness. Neat.
        =#
    end
end

function insphere(pa, pb, pc, pd, pe)
    @inbounds begin
        aex = pa[1] - pe[1]
        bex = pb[1] - pe[1]
        cex = pc[1] - pe[1]
        dex = pd[1] - pe[1]
        aey = pa[2] - pe[2]
        bey = pb[2] - pe[2]
        cey = pc[2] - pe[2]
        dey = pd[2] - pe[2]
        aez = pa[3] - pe[3]
        bez = pb[3] - pe[3]
        cez = pc[3] - pe[3]
        dez = pd[3] - pe[3]

        aexbey = aex * bey
        bexaey = bex * aey
        ab = aexbey - bexaey
        bexcey = bex * cey
        cexbey = cex * bey
        bc = bexcey - cexbey
        cexdey = cex * dey
        dexcey = dex * cey
        cd = cexdey - dexcey
        dexaey = dex * aey
        aexdey = aex * dey
        da = dexaey - aexdey

        aexcey = aex * cey
        cexaey = cex * aey
        ac = aexcey - cexaey
        bexdey = bex * dey
        dexbey = dex * bey
        bd = bexdey - dexbey

        abc = aez * bc - bez * ac + cez * ab
        bcd = bez * cd - cez * bd + dez * bc
        cda = cez * da + dez * ac + aez * cd
        dab = dez * ab + aez * bd + bez * da

        alift = aex * aex + aey * aey + aez * aez
        blift = bex * bex + bey * bey + bez * bez
        clift = cex * cex + cey * cey + cez * cez
        dlift = dex * dex + dey * dey + dez * dez

        det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd)

        aezplus = Absolute(aez)
        bezplus = Absolute(bez)
        cezplus = Absolute(cez)
        dezplus = Absolute(dez)
        aexbeyplus = Absolute(aexbey)
        bexaeyplus = Absolute(bexaey)
        bexceyplus = Absolute(bexcey)
        cexbeyplus = Absolute(cexbey)
        cexdeyplus = Absolute(cexdey)
        dexceyplus = Absolute(dexcey)
        dexaeyplus = Absolute(dexaey)
        aexdeyplus = Absolute(aexdey)
        aexceyplus = Absolute(aexcey)
        cexaeyplus = Absolute(cexaey)
        bexdeyplus = Absolute(bexdey)
        dexbeyplus = Absolute(dexbey)
        permanent = ((cexdeyplus + dexceyplus) * bezplus +
                     (dexbeyplus + bexdeyplus) * cezplus +
                     (bexceyplus + cexbeyplus) * dezplus) *
                    alift +
                    ((dexaeyplus + aexdeyplus) * cezplus +
                     (aexceyplus + cexaeyplus) * dezplus +
                     (cexdeyplus + dexceyplus) * aezplus) *
                    blift +
                    ((aexbeyplus + bexaeyplus) * dezplus +
                     (bexdeyplus + dexbeyplus) * aezplus +
                     (dexaeyplus + aexdeyplus) * bezplus) *
                    clift +
                    ((bexceyplus + cexbeyplus) * aezplus +
                     (cexaeyplus + aexceyplus) * bezplus +
                     (aexbeyplus + bexaeyplus) * cezplus) *
                    dlift
        errbound = isperrboundA(permanent) * permanent
        if (det > errbound) || (-det > errbound)
            return det
        end

        return insphereadapt(pa, pb, pc, pd, pe, permanent)
    end
end