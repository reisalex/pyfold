"""
Module: DELTAG_HD (RNA,II,JJ,KK)

Description: Computes the difference in free energy of an RNA loop due
             to a nucleotide diffusion along the helix i.e. the ii-jj
             base pair will shift to either ii-kk or kk-jj. The energy
             change is calculated using the empirical INN model.

Arguments:

            RNA - Class containing information about the RNA fold and
                the RNA sequence.
             II - Nucleotide position of the 5' most nucleotide.
             JJ - Nucleotide position of the 3' most nucleotide.
             KK - Nucleotide position of the single-stranded nucleotiode
                that either ii/jj in the base-pair ii-jj will swap with.

History:
Version     Date            Comment
--------    -------         --------------------
            09/28/2017      Original Code

Dependencies:

Author(s): Alex Reis
             Copyright (c) 2017 (Please refer to LICENCE)
"""

import math
from rnavar import em,eh,es,eaup,beta

def DELTAG_HD(rna,ii,jj,kk):

    # INTEGER
    # ii,jj,kk
    # i,j,k,l,n,ip,jp
    # nh,ns,mh,hs,js
    # indx,indx2,i1,iloop,kloop
    #note replaced is with hs for Cython code

    # REAL
    # e1,e2,e3,e4,e5,ed
    # x,c,ei,ef

    indx = rna.link[jj]

    i = rna.loop[indx]
    j = rna.ibsp[i]
    n = rna.n

    nh = rna.nhlx[indx]
    ns = rna.nsgl[indx]

    c = 1.750e0 / float(beta)

    #=== Internal Loop iloop = 1 ===#
    #=== External Loop iloop = 0 ===#

    iloop = 1
    kloop = 1

    if ( i > j ):
        iloop = 0
        i = 0
        j = n-1

    indx2 = -1

    if ( rna.link[ii] != -1 ): # ii-jj defines another loop-k, next to loop-indx

        indx2 = rna.link[ii]

        mh = rna.nhlx[indx2]
        ms = rna.nsgl[indx2]

        k = rna.loop[indx2]
        l = rna.ibsp[k]

        if ( k > l ):
            kloop = 0
            # k = 0
            # l = n-1

    else: # ii-jj is part of a helix

        mh = 2
        ms = 0

        # k = ii # not used
        # l = jj

    e4 = 0.0e0

    #=== INITIAL ENERGY ===#

    # First loop

    if ( nh == 1 and iloop == 1 ):

        e1 = EHAIR(rna.iseq,i,j,n)

    elif ( nh == 2 and iloop == 1 ):

        ip = i + 1
        while ( rna.ibsp[ip] == -1 ):
            ip += 1

        jp = rna.ibsp[ip]

        e1 = EBULGE(rna.iseq,i,j,ip,jp,n)

    else:

        e1 = 0.0e0
        e5 = 0.0e0
        e3 = 0.0e0

        hs = rna.iseq[ii]
        js = rna.iseq[jj]

        if ( iloop == 1 ):

            if ( params.MBLmodel == 2 ):

                # calculate average asymmetry of MBL
                asym = 0.0e0
                for i1 in xrange(1,nh+1):
                    asym += float( abs( rna.lns[indx][i1] - rna.lns[indx][i1-1] ) )
                asym /= float(nh)
                asym = min(2.0,asym)

                e1  = params.MBLinit[0]               # a
                e1 += params.MBLinit[1] * float(asym) # b
                e1 += params.MBLinit[2] * float(nh)   # c

                if (nh == 3) and (ns < 2):
                    e1 += params.MBLinit[4] # dG_strain

            elif ( ns <= 6 ):
                e1  = params.MBLinit[0]               # a, (em)
                e1 += params.MBLinit[1] * float(ns)   # b, (es)
                e1 += params.MBLinit[2] * float(nh)   # c, (eh)

            else:
                e1  = params.MBLinit[0]
                e1 += params.MBLinit[1] * 6.0e0
                x   = float(ns) / 6.0e0
                e1 += c * math.log(x)
                e1 += params.MBLinit[2] * float(nh)


        if ( ii > 0 ) and ( rna.ibsp[ii-1] == -1 ):

            e5 = EDANGLE(rna.iseq,ii,jj,ii-1,n)

            if ( ii > 1 ):

                ip = ii - 2
                jp = rna.ibsp[ip]

                if ( jp != -1 ):
                    ed = EDANGLE(rna.iseq,ip,jp,ii-1,n)
                    if ( kk == ii+1 ): e4 += ed
                    e5 = min(e5,ed)

            if  ( ii > 2 and kk == ii-1 ) and
                ( rna.ibsp[ii-2] == -1 ):

                ip = ii - 3
                jp = rna.ibsp[ip]

                if ( jp != -1 ):
                    ed = EDANGLE(rna.iseq,ip,jp,ii-2,n)
                    e1 += ed

        if ( jj < n-1 ) and ( rna.ibsp[jj+1] == -1 ):

            e3 = EDANGLE(rna.iseq,ii,jj,jj+1,n)

            if ( jj < n-2 ):

                ip = jj + 2
                jp = rna.ibsp[ip]

                if ( jp != -1 ):
                    ed = EDANGLE(rna.iseq,ip,jp,jj+1,n)
                    if ( kk == jj-1 ): e4 += ed
                    e3 = min(e3,ed)

            if  ( jj < n-3 and kk == jj+1 ) and
                ( rna.ibsp[jj+2] == -1 ):

                ip = jj + 3
                jp = rna.ibsp[ip]

                if ( jp != -1 ):
                    ed = EDANGLE(rna.iseq,ip,jp,jj+2,n)
                    e1 += ed

        e1 = e1 + e3 + e5 + params.dG_AUP[hs][js]

    # Second loop

    if ( mh == 1 and kloop == 1 ):

        e2 = EHAIR(rna.iseq,ii,jj,n)

    elif ( mh == 2 and kloop == 1 ):

        if ( ms == 0 ):

            e2 = ESTACK(rna.iseq,ii,jj,ii+1,jj-1,n)

        else:

            ip = ii + 1
            while ( rna.ibsp[ip] == -1 ):
                ip += 1

            jp = rna.ibsp[ip]

            e2 = EBULGE(rna.iseq,ii,jj,ip,jp,n)

    else:

        e2 = 0.0e0
        e5 = 0.0e0
        e3 = 0.0e0

        hs = rna.iseq[ii]
        js = rna.iseq[jj]

        if ( kloop == 1 ):

            if ( params.MBLmodel == 2 ):

                # calculate average asymmetry of MBL
                asym = 0.0e0
                for i1 in xrange(1,mh+1):
                    asym += float( abs( rna.lns[indx2][i1] - rna.lns[indx2][i1-1] ) )
                asym /= float(mh)
                asym = min(2.0,asym)

                e2  = params.MBLinit[0]               # a
                e2 += params.MBLinit[1] * float(asym) # b
                e2 += params.MBLinit[2] * float(mh)   # c

                if (mh == 3) and (ms < 2):
                    e2 += params.MBLinit[4] # dG_strain

            elif ( ms <= 6 ):
                e2  = params.MBLinit[0]               # a, (em)
                e2 += params.MBLinit[1] * float(ms)   # b, (es)
                e2 += params.MBLinit[2] * float(mh)   # c, (eh)

            else:
                e2  = params.MBLinit[0]
                e2 += params.MBLinit[1] * 6.0e0
                x   = float(ms) / 6.0e0
                e2 += c * math.log(x)
                e2 += params.MBLinit[2] * float(mh)

        if ( jj > 0 ) and ( rna.ibsp[jj-1] == -1 ):

            e5 = EDANGLE(rna.iseq,ii,jj,jj-1,n)

            if ( jj > 1 ):

                ip = jj - 2
                jp = rna.ibsp[ip]

                if ( jp != -1 ):
                    ed = EDANGLE(rna.iseq,ip,jp,jj-1,n)
                    if ( kk == jj+1 ): e4 += ed
                    e5 = min(e5,ed)

            if  ( jj > 2 and kk == jj-1 ) and
                ( rna.iseq[jj-2] == -1 ):

                ip = jj - 3
                jp = rna.ibsp[ip]

                if ( jp != -1 ):
                    ed = EDANGLE(rna.iseq,ip,jp,jj-2,n)
                    e2 += ed

        if ( ii < n-1 ) and ( rna.ibsp[ii+1] == -1 ):

            e3 = EDANGLE(rna.iseq,ii,jj,ii+1,n)

            if ( ii < n-2 ):

                ip = ii + 2
                jp = rna.ibsp[ip]

                if ( jp != -1 ):
                    ed = EDANGLE(rna.iseq,ip,jp,ii+1,n)
                    if ( kk == ii-1 ): e4 += ed
                    e3 = min(e3,ed)

            if  ( ii < n-3 and kk == ii+1 ) and
                ( rna.ibsp[ii+2] == -1 ):

                ip = ii + 3
                jp = rna.ibsp[ip]

                if ( jp != -1 ):
                    ed = EDANGLE(rna.iseq,ip,jp,ii+2,n)
                    e2 += ed

        e2 = e2 + e3 + e5 + params.dG_AUP[hs][js]

    ei = e1 + e2

    #=== FINAL ENERGY ===#

    # lns for loop i-j
    lns1 = rna.lns[indx][:]
    hi1  = rna.htrack[indx][ii]

    if kk == ii-1:
        ns -= 1
        ms += 1
        lns1[hi1-1] -= 1

        if indx2 != -1:
            lns2 = rna.lns[indx2][:]
            lns2[1] += 1

    elif kk == jj+1:
        ns -= 1
        ms += 1
        lns1[hi1] -= 1

        if indx2 != -1:
            lns2 = rna.lns[indx2][:]
            lns2[mh] += 1

    elif kk == ii+1:
        ns += 1
        ms -= 1
        lns1[hi1-1] += 1

        if indx2 != -1:
            lns2 = rna.lns[indx2][:]
            lns2[1] -= 1

    else: # kk == jj-1
        ns += 1
        ms -= 1
        lns1[hi1] += 1

        if indx2 != -1:
            lns2 = rna.lns[indx2][:]
            lns2[mh] -= 1

    # First loop

    if ( nh == 1 and iloop == 1 ):

        if ( kk == i-1 or kk == i+1 ):
            e1 = EHAIR(rna.iseq,kk,j,n)

        if ( kk == j-1 or kk == j+1 ):
            e1 = EHAIR(rna.iseq,i,kk,n)

    elif ( nh == 2 and iloop == 1 ):

        if ( ns == 0 ):

            if ( kk == ii-1 ):
                e1 = ESTACK(rna.iseq,kk-1,jj+1,kk,jj,n)
            else:
                e1 = ESTACK(rna.iseq,ii-1,kk+1,ii,kk,n)

        else:

            ip = i + 1
            while ( rna.ibsp[ip] == -1 ):
                ip += 1

            jp = rna.ibsp[ip]

            if ( ii == j ):

                if ( kk == ii+1 or kk == ii-1 ):
                    e1 = EBULGE(rna.iseq,jj,kk,ip,jp,n)
                else:
                    e1 = EBULGE(rna.iseq,kk,ii,ip,jp,n)

            else:

                if ( kk == ii+1 or kk == ii-1 ):
                    e1 = EBULGE(rna.iseq,i,j,kk,jj,n)
                else:
                    e1 = EBULGE(rna.iseq,i,j,ii,kk,n)

    else:

        e1 = 0.0e0
        e5 = 0.0e0
        e3 = 0.0e0

        if ( iloop == 1 ):

            if ( params.MBLmodel == 2 ):

                # calculate average asymmetry of MBL
                asym = 0.0e0
                for i1 in xrange(1,nh+1):
                    asym += float( abs( lns1[i1] - lns1[i1-1] ) )
                asym /= float(nh)
                asym = min(2.0,asym)

                e1  = params.MBLinit[0]               # a
                e1 += params.MBLinit[1] * float(asym) # b
                e1 += params.MBLinit[2] * float(nh)   # c

                if (nh == 3) and (ns < 2):
                    e1 += params.MBLinit[4] # dG_strain

            elif ( ns <= 6 ):
                e1  = params.MBLinit[0]               # a, (em)
                e1 += params.MBLinit[1] * float(ns)   # b, (es)
                e1 += params.MBLinit[2] * float(nh)   # c, (eh)

            else:
                e1  = params.MBLinit[0]
                e1 += params.MBLinit[1] * 6.0e0
                x   = float(ns) / 6.0e0
                e1 += c * math.log(x)
                e1 += params.MBLinit[2] * float(nh)

        if ( kk == jj+1 or kk == jj-1 ):

            hs = rna.iseq[ii]
            js = rna.iseq[kk]

            if ( ii > 0 ) and ( rna.ibsp[ii-1] == -1 ):

                e5 = EDANGLE(rna.iseq,ii,kk,ii-1,n)

                if ( ii > 1 ):

                    ip = ii - 2
                    jp = rna.ibsp[ip]

                    if ( jp != -1 ):
                        ed = EDANGLE(rna.iseq,ip,jp,ii-1,n)
                        e5 = min(e5,ed)

            if ( kk < n-1 ) and ( rna.ibsp[kk+1] == -1 or kk == jj-1 ):

                e3 = EDANGLE(rna.iseq,ii,kk,kk+1,n)

                if ( kk < n-2 ):

                    ip = kk + 2
                    jp = rna.ibsp[ip]

                    if ( jp != -1 ):
                        ed = EDANGLE(rna.iseq,ip,jp,kk+1,n)
                        e3 = min(e3,ed)
        
        if ( kk == ii-1 or kk == ii+1 ):

            hs = rna.iseq[kk]
            js = rna.iseq[jj]

            if ( kk > 0 ) and (rna.ibsp[kk-1] == -1 or kk == ii+1 ):

                e5 = EDANGLE(rna.iseq,kk,jj,kk-1,n)

                if ( kk > 1 ):

                    ip = kk - 2
                    jp = rna.ibsp[ip]

                    if ( jp != -1 ):
                        ed = EDANGLE(rna.iseq,ip,jp,kk-1,n)
                        e5 = min(e5,ed)

            if ( jj < n-1 ) and ( rna.ibsp[jj+1] == -1 ):

                e3 = EDANGLE(rna.iseq,kk,jj,jj+1,n)

                if ( jj < n-2 ):

                    ip = jj + 2
                    jp = rna.ibsp[ip]

                    if ( jp != -1 ):
                        ed = EDANGLE(rna.iseq,ip,jp,jj+1,n,ed)
                        e3 = min(e3,ed)

        e1 = e1 + e3 + e5 + params.dG_AUP[hs][js]

    # Second loop

    if ( mh == 1 and kloop == 1 ):

        if ( kk == ii-1 or kk == ii+1 ):
            e2 = EHAIR(rna.iseq,kk,jj,n)

        if ( kk == jj-1 or kk == jj+1 ):
            e2 = EHAIR(rna.iseq,ii,kk,n)

    elif ( mh == 2 and kloop == 1 ):

        if ( ms == 0 ):

            if ( kk == ii+1 ):
                e2 = ESTACK(rna.iseq,kk,jj,kk+1,jj-1,n)
            else:
                e2 = ESTACK(rna.iseq,ii,kk,ii+1,kk-1,n)
        
        else:

            ip = ii + 1
            while ( rna.ibsp[ip] == -1 ):
                ip += 1

            jp = rna.ibsp[ip]

            if ( kk == ii+1 or kk == ii-1 ):
                e2 = EBUGLE(rna.iseq,kk,jj,ip,jp,n)
            else:
                e2 = EBUGLE(rna.iseq,ii,kk,ip,jp,n)

    else:

        e2 = 0.0e0
        e5 = 0.0e0
        e3 = 0.0e0

        if ( kloop == 1 ):

            if ( params.MBLmodel == 2 ):

                # calculate average asymmetry of MBL
                asym = 0.0e0
                for i1 in xrange(1,mh+1):
                    asym += float( abs( lns2[i1] - lns2[i1-1] ) )
                asym /= float(mh)
                asym = min(2.0,asym)

                e2  = params.MBLinit[0]               # a
                e2 += params.MBLinit[1] * float(asym) # b
                e2 += params.MBLinit[2] * float(mh)   # c

                if (mh == 3) and (ms < 2):
                    e2 += params.MBLinit[4] # dG_strain

            elif ( ms <= 6 ):
                e2  = params.MBLinit[0]               # a, (em)
                e2 += params.MBLinit[1] * float(ms)   # b, (es)
                e2 += params.MBLinit[2] * float(mh)   # c, (eh)

            else:
                e2  = params.MBLinit[0]
                e2 += params.MBLinit[1] * 6.0e0
                x   = float(ms) / 6.0e0
                e2 += c * math.log(x)
                e2 += params.MBLinit[2] * float(mh)

        if ( kk == jj+1 or kk == jj-1 ):

            hs = rna.iseq[ii]
            js = rna.iseq[kk]

            if ( kk > 0 ) and ( rna.ibsp[kk-1] == -1 or kk == jj+1 ):

                e5 = EDANGLE(rna.iseq,ii,kk,kk-1,n)

                if ( kk > 1 ):

                    ip = kk - 2
                    jp = rna.ibsp[ip]

                    if ( jp != -1 ):
                        ed = EDANGLE(rna.iseq,ip,jp,kk-1,n)
                        e5 = min(e5,ed)

            if ( ii < n-1 ) and ( rna.ibsp[ii+1] == -1 ):

                e3 = EDANGLE(rna.iseq,ii,kk,ii+1,n)

                if ( ii < n-2 ):

                    ip = ii + 2
                    jp = rna.ibsp[ip]

                    if ( jp != -1 ):
                        ed = EDANGLE(rna.iseq,ip,jp,ii+1,n)
                        e3 = min(e3,ed)

        if ( kk == ii-1 or kk == ii+1 ):

            hs = rna.iseq[kk]
            js = rna.iseq[jj]

            if ( jj > 0 ) and ( rna.ibsp[jj-1] == -1 ):

                e5 = EDANGLE(rna.iseq,kk,jj,jj-1,n)

                if ( jj > 1 ):

                    ip = jj - 2
                    jp = rna.ibsp[ip]

                    if ( jp != -1 ):
                        ed = EDANGLE(rna.iseq,ip,jp,jj-1,n)
                        e5 = min(e5,ed)

            if ( kk < n-1 ) and ( rna.ibsp[kk+1] == -1 or kk == ii-1 ):

                e3 = EDANGLE(rna.iseq,kk,jj,kk+1,n)

                if ( kk < n-1 ):

                    ip = kk + 2
                    jp = rna.ibsp[ip]

                    if ( jp != -1 ):
                        ed = EDANGLE(rna.iseq,ip,jp,kk+1,n)
                        e3 = min(e3,ed)

        e2 = e2 + e3 + e5 + params.dG_AUP[hs][js]

    #endif

    ef = e1 + e2 + e4

    dg = float(ef) - float(ei)

    return dg





