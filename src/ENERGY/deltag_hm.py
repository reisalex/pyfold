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
    # is,js,nh,ns,mh,ms

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

    if ( i > j ):
        iloop = 0
        i = 1
        j = n
    #endif

    #=== Final Loop Size ===#

    mh = nh
    ms = ns

    ip = rna.ibsp[ii-1]
    jp = rna.ibsp[jj+1]

    if ( ip != 0 and jp != 0 ): ms += 2

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

        if ( ns <= 6 ):
            ei = em + es * float(ns) + eh * float(nh)
        else:
            x = float(ns) / 6.0e0
            ei = em + es * 6.0e0 + eh * float(nh)
            ei = ei + c * math.log(x)
        #endif
    #endif

    ei += eaup[hs][js]

    #=== 5' Side ===#

    if ( rna.ibsp[ii-1] != 0 ):

        jp = ii - 1
        ip = rna.ibsp[jp]

        hs = rna.iseq[ip]
        hs = rna.iseq[jp]

        e1 = ESTACK(rna.iseq,ip,jp,ip+1,jp-1,n)

        if ( ip > 1 ) and ( rna.ibsp[ip-1] == 0 ):

            e5 = EDANGLE(rna.iseq,ip,jp,ip-1,n)

            if ( ip > 2 ):

                kp = ip - 2
                lp = rna.ibsp[kp]

                if ( lp != 0 ):
                    ed = EDANGLE(rna.iseq,kp,lp,ip-1,n)
                    if ( lp != jj+1 ): e4 += ed
                    e5 = min(e5,ed)
                #endif
            #endif
        #endif

        e1 = e1 + e5 + eaup[hs][js]

        hs = rna.iseq[ip+1]
        js = rna.iseq[jp-1]

        e5 = EDANGLE(rna.iseq,ip+1,jp-1,ip,n)

        if ( ip > 1 ) and ( rna.ibsp[ip-1] != 0 ):

            kp = ip - 1
            lp = rna.ibsp[kp]

            if ( lp != jj+1 ):
                ed = EDANGLE(rna.iseq,kp,lp,ip,n)
                e5 = MIN(e5,ed)
            #endif
        #endif

        e4 = e4 + e5 + eaup[hs][js]

    else:

        e5 = EDANGLE(rna.iseq,ii,jj,ii-1,n)

        if ( ii > 2 ):

            ip = ii - 2
            jp = rna.ibsp[ip]

            if ( jp != 0 ):

                ed = EDANGLE(rna.iseq,ip,jp,ii-1,n)
                e5 = min(e5,ed)

            elif ( ii > 3 ):

                kp = ii - 3
                lp = rna.ibsp[kp]

                if ( lp != 0 ):
                    ed = EDANGLE(rna.iseq,kp,lp,ii-2,n)
                    e5 += ed
                #endif
            #endif
        #endif

        e1 = e5

    #endif

    #=== 3' Side ===#

    if ( rna.ibsp[jj+1] != 0 ):

        ip = jj + 1
        jp = rna.ibsp[ip]

        hs = rna.iseq[ip]
        js = rna.iseq[jp]

        e2 = ESTACK(rna.iseq,ip,jp,ip+1,jp-1,n)

        if ( jp < n ) and ( rna.ibsp[jp+1] == 0 ):

            e3 = EDANGLE(rna.iseq,ip,jp,jp+1,n)

            if ( jp < n-1 ):

                kp = jp + 2
                lp = rna.ibsp[kp]

                if ( lp != 0 ):

                    e3 = EDANGLE(rna.iseq,kp,lp,jp+1,n)

                    if ( lp != ii-1 ):
                        e4 += ed
                        e3 = min(e3,ed)
                    else:
                        e3 = 0.0e0
                    #endif
                #endif
            #endif
        #endif

        e2 = e2 + e3 + eaup[hs][js]

        hs = rna.iseq[ip+1]
        js = rna.iseq[jp-1]

        e3 = EDANGLE(rna.iseq,ip+1,jp-1,jp,n)

        if ( jp < n ) and ( rna.ibsp[jp+1] != 0 ):

            kp = jp + 1
            lp = rna.ibsp[kp]

            if ( lp != ii-1 ):
                ed = EDANGLE(rna.iseq,kp,lp,jp,n)
                e3 = min(e3,ed)
            #endif
        #endif

        e4 = e4 + e3 + eaup[hs][js]

    else:

        e3 = EDANGLE(rna.iseq,ii,jj,jj+1,n)

        if ( jj < n-1 ):

            ip = jj + 2
            jp = rna.ibsp[ip]

            if ( jp != 0 ):

                ed = EDANGLE(rna.iseq,ip,jp,jj+1,n)
                e3 = min(e3,ed)

            elif ( jj < n-2 ):

                kp = jj + 3
                lp = rna.ibsp[kp]

                if ( lp != 0 ):
                    ed = EDANGLE(rna.iseq,kp,lp,jj+2,n)
                    e3 += ed
                #endif
            #endif
        #endif

        e2 = e3

    #endif

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

        if ( ms <= 6 ):
            ef = em + es * float(ms) + eh * float(mh)
        else:
            x = float(ms) / 6.0e0
            ef = em + es * 6.0e0 + eh * float(mh)
            ef = ef + c * math.log(x)
        #endif
    #endif

    ed = ESTACK(rna.iseq,ip,jp,ii,jj,n)

    ef = ef + ed + eaup[hs][js]

    if ( ip > 1 ) and ( rna.ibsp[ip-1] == 0 ):

        e5 = EDANGLE(rna.iseq,ip,jp,ip-1,n)

        if ( ip > 2 ):

            kp = ip - 2
            lp = rna.ibsp[kp]

            if ( lp != 0 ):
                ed = EDANGLE(rna.iseq,kp,lp,ip-1,n)
                e5 = min(e5,ed)
            #endif
        #endif
    #endif

    if ( jp < n ) and ( rna.ibsp[jp+1] == 0 ):

        e3 = EDANGLE(rna.iseq,ip,jp,jp+1,n)

        if ( jp < n-1 ):

            kp = jp + 2
            lp = rna.ibsp[kp]

            if ( lp != 0 ):
                ed = EDANGLE(rna.iseq,kp,lp,jp+1,n)
                e3 = min(e3,ed)
            #endif
        #endif
    #endif

    ef = ef + e5 + e3 + e4

    dg = float(ef) - float(ei)

    return dg
