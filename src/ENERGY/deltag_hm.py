"""
Module: DELTAG_HM (RNA,II,JJ)

Description: Computes the difference in free energy of an RNA loop due to
             the swapping of base-pairs between two helices using the
             empirical INN model. The swapping of base-pairs occurs as a
             result of the extension of the helix below the ii-jj base-pair.

Arguments:

            R - Class containing information about the RNA fold and
                the RNA sequence.
             II - Nucleotide position of the 5' most nucleotide.
             JJ - Nucleotide position of the 3' most nucleotide.

History:
Version     Date            Comment
--------    -------         --------------------
            09/28/2017      Original Code

Dependencies:

Author(s): Alex Reis
             Copyright (c) 2017 (Please refer to LICENCE)
"""

def DELTAG_HM(rna,ii,jj):

    # INTEGER
    # ii,jj
    # i,j,k,n,ip,jp,kp,lp
    # is,js,nh,ns,mh,ms,hi,lns
    # b1,b2,i1

    # REAL
    # e1,e2,e3,e4,e5,ed
    # x,c,ei,ef

    indx = rna.link[jj]

    i = rna.loop[indx]
    j = rna.ibsp[i]
    n = rna.n

    nh = rna.nhlx[indx]
    ns = rna.nsgl[indx]

    hi = rna.htrack[indx][ii]
    lns = rna.lns[indx][:]

    c = 1.750e0 / float(beta)

    #=== Internal Loop iloop = 1 ===#
    #=== External Loop iloop = 0 ===#

    iloop = 1

    if ( i > j ):
        iloop = 0
        i = 0
        j = n-1

    #=== Final Loop Size ===#

    mh = nh
    ms = ns

    ip = rna.ibsp[ii-1]
    jp = rna.ibsp[jj+1]

    # are there adjacent helices?
    b1 = ( ip != -1 )
    b2 = ( jp != -1 )

    if ( b1 and b2 ):
        # free up two base pairs in loop indx
        ms += 2
        lns[hi-2] += 1
        lns[hi+1] += 1
    else:
        if not b1: # remove 1 nt from lns
            lns[hi-1] -= 1
        else: # by helix morphing definition
            lns[hi] -= 1

    e4 = 0.0e0

    #=== Initial Energy ===#

    ei = 0.0e0
    e1 = 0.0e0
    e2 = 0.0e0
    e3 = 0.0e0
    e5 = 0.0e0

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

            ei  = params.MBLinit[0]               # a
            ei += params.MBLinit[1] * float(asym) # b
            ei += params.MBLinit[2] * float(nh)   # c

            if (nh == 3) and (ns < 2):
                ei += params.MBLinit[4] # dG_strain

        elif ( ns <= 6 ):
            ei  = params.MBLinit[0]               # a, (em)
            ei += params.MBLinit[1] * float(ns)   # b, (es)
            ei += params.MBLinit[2] * float(nh)   # c, (eh)

        else:
            ei  = params.MBLinit[0]
            ei += params.MBLinit[1] * 6.0e0
            x   = float(ns) / 6.0e0
            ei += c * math.log(x)
            ei += params.MBLinit[2] * float(nh)

    ei += params.dG_AUP[hs][js]

    #=== 5' Side ===#

    if ( rna.ibsp[ii-1] != -1 ):

        jp = ii - 1
        ip = rna.ibsp[jp]

        hs = rna.iseq[ip]
        js = rna.iseq[jp]

        e1 = ESTACK(rna.iseq,ip,jp,ip+1,jp-1,n)

        if ( ip > 0 ) and ( rna.ibsp[ip-1] == -1 ):

            e5 = EDANGLE(rna.iseq,ip,jp,ip-1,n)

            if ( ip > 1 ):

                kp = ip - 2
                lp = rna.ibsp[kp]

                if ( lp != -1 ):
                    ed = EDANGLE(rna.iseq,kp,lp,ip-1,n)
                    if ( lp != jj+1 ): e4 += ed
                    e5 = min(e5,ed)

        e1 = e1 + e5 + params.dG_AUP[hs][js]

        hs = rna.iseq[ip+1]
        js = rna.iseq[jp-1]

        e5 = EDANGLE(rna.iseq,ip+1,jp-1,ip,n)

        if ( ip > 0 ) and ( rna.ibsp[ip-1] != -1 ):

            kp = ip - 1
            lp = rna.ibsp[kp]

            if ( lp != jj+1 ):
                ed = EDANGLE(rna.iseq,kp,lp,ip,n)
                e5 = MIN(e5,ed)

        e4 = e4 + e5 + params.dG_AUP[hs][js]

    else:

        e5 = EDANGLE(rna.iseq,ii,jj,ii-1,n)

        if ( ii > 1 ):

            ip = ii - 2
            jp = rna.ibsp[ip]

            if ( jp != -1 ):

                ed = EDANGLE(rna.iseq,ip,jp,ii-1,n)
                e5 = min(e5,ed)

            elif ( ii > 2 ):

                kp = ii - 3
                lp = rna.ibsp[kp]

                if ( lp != -1 ):
                    ed = EDANGLE(rna.iseq,kp,lp,ii-2,n)
                    e5 += ed

        e1 = e5

    #=== 3' Side ===#

    if ( rna.ibsp[jj+1] != -1 ):

        ip = jj + 1
        jp = rna.ibsp[ip]

        hs = rna.iseq[ip]
        js = rna.iseq[jp]

        e2 = ESTACK(rna.iseq,ip,jp,ip+1,jp-1,n)

        if ( jp < n-1 ) and ( rna.ibsp[jp+1] == -1 ):

            e3 = EDANGLE(rna.iseq,ip,jp,jp+1,n)

            if ( jp < n-2 ):

                kp = jp + 2
                lp = rna.ibsp[kp]

                if ( lp != -1 ):

                    ed = EDANGLE(rna.iseq,kp,lp,jp+1,n)

                    if ( lp != ii-1 ):
                        e4 += ed
                        e3 = min(e3,ed)
                    else:
                        e3 = 0.0e0

        e2 = e2 + e3 + params.dG_AUP[hs][js]

        hs = rna.iseq[ip+1]
        js = rna.iseq[jp-1]

        e3 = EDANGLE(rna.iseq,ip+1,jp-1,jp,n)

        if ( jp < n-1 ) and ( rna.ibsp[jp+1] != -1 ):

            kp = jp + 1
            lp = rna.ibsp[kp]

            if ( lp != ii-1 ):
                ed = EDANGLE(rna.iseq,kp,lp,jp,n)
                e3 = min(e3,ed)

        e4 = e4 + e3 + params.dG_AUP[hs][js]

    else:

        e3 = EDANGLE(rna.iseq,ii,jj,jj+1,n)

        if ( jj < n-2 ):

            ip = jj + 2
            jp = rna.ibsp[ip]

            if ( jp != -1 ):

                ed = EDANGLE(rna.iseq,ip,jp,jj+1,n)
                e3 = min(e3,ed)

            elif ( jj < n-3 ):

                kp = jj + 3
                lp = rna.ibsp[kp]

                if ( lp != -1 ):
                    ed = EDANGLE(rna.iseq,kp,lp,jj+2,n)
                    e3 += ed

        e2 = e3

    ei = ei + e1 + e2

    #=== FINAL ENERGY ===#

    ef = 0.0e0
    e5 = 0.0e0
    e3 = 0.0e0

    ip = ii - 1
    jp = jj + 1

    hs = rna.iseq[ip]
    js = rna.iseq[jp]

    if ( iloop == 1 ):

        if ( params.MBLmodel == 2 ):

            # calculate average asymmetry of MBL
            asym = 0.0e0
            for i1 in xrange(1,mh+1):
                asym += float( abs( lns[i1] - lns[i1-1] ) )
            asym /= float(mh)
            asym = min(2.0,asym)

            ef  = params.MBLinit[0]               # a
            ef += params.MBLinit[1] * float(asym) # b
            ef += params.MBLinit[2] * float(mh)   # c

            if (mh == 3) and (ms < 2):
                ef += params.MBLinit[4] # dG_strain

        elif ( ms <= 6 ):
            ef  = params.MBLinit[0]               # a, (em)
            ef += params.MBLinit[1] * float(ms)   # b, (es)
            ef += params.MBLinit[2] * float(mh)   # c, (eh)

        else:
            ef  = params.MBLinit[0]
            ef += params.MBLinit[1] * 6.0e0
            x   = float(ms) / 6.0e0
            ef += c * math.log(x)
            ef += params.MBLinit[2] * float(mh)

    ed = ESTACK(rna.iseq,ip,jp,ii,jj,n)

    ef = ef + ed + params.dG_AUP[hs][js]

    if ( ip > 0 ) and ( rna.ibsp[ip-1] == -1 ):

        e5 = EDANGLE(rna.iseq,ip,jp,ip-1,n)

        if ( ip > 1 ):

            kp = ip - 2
            lp = rna.ibsp[kp]

            if ( lp != -1 ):
                ed = EDANGLE(rna.iseq,kp,lp,ip-1,n)
                e5 = min(e5,ed)

    if ( jp < n-1 ) and ( rna.ibsp[jp+1] == -1 ):

        e3 = EDANGLE(rna.iseq,ip,jp,jp+1,n)

        if ( jp < n-2 ):

            kp = jp + 2
            lp = rna.ibsp[kp]

            if ( lp != -1 ):
                ed = EDANGLE(rna.iseq,kp,lp,jp+1,n)
                e3 = min(e3,ed)

    ef = ef + e5 + e3 + e4

    dg = float(ef) - float(ei)

    return dg
