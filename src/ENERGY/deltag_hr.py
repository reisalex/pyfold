"""
Function: DELTAG_HR (RNA,II,JJ)

Description: Computes the difference in free energy of an RNA loop due to
             a deletion of the ii-jj base pair using the empirical INN model.

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

import math
from rnavar import em,eh,es,eaup,beta

def DELTAG_HR(rna,ii,jj):

    # VARIABLES
    
    # INTEGER
    # ii,jj
    # i,j,k,n,ip,jp,hs,js
    # nh,ns,mh,ms,lh,ls
    # iloop, indx, lnsf, hi, hjmp, lnsf, i1

    # DOUBLE PRECISION
    # dg

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
        i = 0
        j = n-1

    if ( rna.ibsp[ii+1] != jj-1 ):

        k  = rna.link[ii]

        mh = rna.nhlx[k]
        ms = rna.nsgl[k]

    else:

        mh = 2
        ms = 0

    #=== Final Loop Size ===#

    lnsf = rna.lns[indx][:]
    hi = rna.htrack[indx][ii]

    if ( mh == 1 ):

        lnsf[hi-1] += 2 + ms + lnsf[hi]

        for i1 in xrange(hi,nh):
            lnsf[i1] = lnsf[i1+1]

    elif ( mh == 2 and ms > 0 ):

        lnsf[hi-1] += rna.lns[k][1] + 1
        lnsf[hi]   += rna.lns[k][0] + 1

    elif ( mh > 2 ):

        hjmp = 0
        
        # add helix parts from loop-indx
        for i1 in xrange(1,nh+1):
            if i1 == hi: hjmp = mh - 2
            lnsf[i1+hjmp] = lnsf[i1]

        # add helix parts from loop-k
        lnsf[hi-1]    += 1 + rna.lns[k][1]
        lnsf[hi+hjmp] += 1 + rna.lns[k][0]

        for i1 in xrange(2,mh):
            lnsf[hi+i1-hjmp] = rna.lnsf[k][i1]

    else: # mh = 2 and ms = 0

        lnsf[hi-1] += 1
        lnsf[hi]   += 1

    lh = nh + mh - 2
    ls = ns + ms + 2

    e4 = 0.0e0

    #=== Initial Energy ===#

    if ( nh == 1 and iloop == 1 ):

        e1 = EHAIR(rna.iseq,i,j,n)

    elif ( nh == 2 and iloop == 1 ):

        ip = i + 1
        while ( rna.ibsp[ip] == 0 ):
            ip += 1
        #endwhile

        jp = rna.ibsp[ip]

        e1 = EBULGE(rna.iseq,i,j,ip,jp,n)

        if ( mh > 2 ):

            hs = rna.iseq[i]
            js = rna.iseq[j]

            if ( rna.ibsp[i+1] == 0 ):
                ed = EDANGLE(rna.iseq,i,j,i+1,n)
                e4 += ed
            #endif

            if ( rna.ibsp[j-1] == 0 ):
                ed = EDANGLE(rna.iseq,i,j,j-1,n)
                e4 += ed
            #endif

            e4 += eaup[hs][js]

        #endif

    else:

        e1 = 0.0e0
        e5 = 0.0e0
        e3 = 0.0e0

        hs = rna.iseq[ii]
        js = rna.iseq[jj]

        if ( iloop == 1 ):

            if ( ns <= 6 ):
                e1 = em + es * float(ns) + eh * float(nh)
            else:
                x = float(ns) / 6.0e0
                e1 = em + es * 6.0e0 + eh * float(nh)
                e1 = e1 + c * math.log(x)
            #endif
        #endif

        if ( ii > 1 ) and ( rna.ibsp[ii-1] == 0 ):

            e5 = EDANGLE(rna.iseq,ii,jj,ii-1,n)

            if ( ii > 2 ):

                ip = ii - 2
                jp = rna.ibsp[ip]

                if ( jp != 0 ):
                    ed = EDANGLE(rna.iseq,ip,jp,ii-1,n)
                    e5 = min(e5,ed)
                    e4 += ed
                #endif
            #endif
        #endif

        if ( jj < n ) and ( rna.ibsp[jj+1] == 0 ):

            e3 = EDANGLE(rna.iseq,ii,jj,jj+1,n)

            if ( jj < n-1 ):

                ip = jj + 2
                jp = rna.ibsp[ip]

                if ( jp != 0 ):
                    ed = EDANGLE(rna.iseq,ip,jp,jj+1,n)
                    e3 = min(e3,ed)
                    e4 += ed
                #endif
            #endif
        #endif

        e1 = e1 + e3 + e5 + eaup[hs][js]

    #endif

    if ( mh == 1 ):

        e2 = EHAIR(rna.iseq,ii,jj,n)

    elif ( mh == 2 ):

        if ( ms == 0 ):

            e2 = ESTACK(rna.iseq,ii,jj,ii+1,jj-1,n)

        else:

            ip = ii + 1
            while ( rna.ibsp[ip] == 0 ):
                ip += 1
            #endwhile

            jp = rna.ibsp[ip]

            e2 = EBULGE(rna.iseq,ii,jj,ip,jp,n)

            if ( nh > 2 or iloop == 0 ):

                hs = rna.iseq[ip]
                js = rna.iseq[jp]

                if ( rna.ibsp[ip-1] == 0 ):
                    ed = EDANGLE(rna.iseq,ip,jp,ip-1,n)
                    e4 += ed
                #endif

                if ( rna.ibsp[jp+1] == 0 ):
                    ed = EDANGLE(rna.iseq,ip,jp,jp+1,n)
                    e4 += ed
                #endif

                e4 += eaup[hs][js]

            #endif
        #endif

    else:

        e2 = 0.0e0
        e5 = 0.0e0
        e3 = 0.0e0

        hs = rna.iseq[ii]
        js = rna.iseq[jj]

        if ( ms <=6 ):
            e2 = em + es * float(ms) + eh * float(mh)
        else:
            x = float(ms) / 6.0e0
            e2 = em + es * 6.0e0 + eh * float(mh)
            e2 = e2 + c * math.log(x)
        #endif

        if ( rna.ibsp[jj-1] == 0 ):

            e5 = EDANGLE(rna.iseq,ii,jj,jj-1,n)

            ip = jj - 2
            jp = rna.ibsp[ip]

            if ( jp != 0 ):
                ed = EDANGLE(rna.iseq,ip,jp,jj-1,n)
                e5 = min(e5,ed)
                e4 += ed
            #endif
        #endif

        if ( rna.ibsp[ii+1] == 0 ):

            e3 = EDANGLE(rna.iseq,ii,jj,ii+1,n)

            ip = ii + 2
            jp = rna.ibsp[ip]

            if ( jp != 0 ):
                ed = EDANGLE(rna.iseq,ip,jp,ii+1,n)
                e3 = min(e3,ed)
                e4 += ed
            #endif
        #endif

        e2 = e2 + e3 + e5 + eaup[hs][js]

    #endif

    ei = e1 + e2

    #=== FINAL ENERGY ===#

    if ( lh == 1 and iloop == 1 ):

        if ( j == ii ):
            ef = EHAIR(rna.iseq,i-1,j+1,n)
        else:
            ef = EHAIR(rna.iseq,i,j,n)
        #endif

    elif ( lh == 2 and iloop == 1 ):

        ip = ii + 1
        if ( mh == 1 ): ip = i + 1
        if ( j == ii ): ip = i + 1

        while ( rna.ibsp[ip] == 0 ):
            ip += 1
        #endwhile

        jp = rna.ibsp[ip]

        if ( nh > 2 ):

            e3 = 0.0e0
            e5 = 0.0e0

            if ( ip == ii ):

                ip = jj + 1
                while ( rna.ibsp[ip] == 0 ):
                    ip += 1
                #endwhile

                jp = rna.ibsp[ip]

            #endif

            hs = rna.iseq[i]
            js = rna.iseq[j]

            if ( i+2 != ii ) and (rna.ibsp[i+1] == 0 ):
                e3 = EDANGLE(rna.iseq,i,j,i+1,n)
            #endif

            if ( j-2 != jj ) and ( rna.ibsp[j-1] == 0 ):
                e5 = EDANGLE(rna.iseq,i,j,j-1,n)
            #endif

            e1 += eaup[hs][js]

            hs = rna.iseq[ip]
            js = rna.iseq[jp]

            if ( ip-2 != jj ) and ( rna.ibsp[ip-1] == 0 ):

                ed = EDANGLE(rna.iseq,ip,jp,ip-1,n)

                if ( ip-2 == i ):
                    e3 = min(e3,ed)
                else:
                    e3 += ed
                #endif

            #endif

            if ( jp+2 != ii ) and ( rna.ibsp[jp+1] == 0 ):

                ed = EDANGLE(rna.iseq,ip,jp,jp+1,n)

                if ( jp+2 == j ):
                    e5 = min(e5,ed)
                else:
                    e5 += ed
                #endif
            #endif

            e1 += eaup[hs][js]
            e1 = e1 + e3 + e5

            ei = e1 + e2 # IS THIS A BUG? SHOULD "ei" be "ef"?

        #endif

        if ( j == ii ):
            ef = EBULGE(rna.iseq,i-1,j+1,ip,jp,n)
        else:
            ef = EBULGE(rna.iseq,i,j,ip,jp,n)
        #endif

    else:

        ef = 0.0e0
        e5 = 0.0e0
        e3 = 0.0e0

        if ( iloop == 1 ):

            if ( ls <= 6 ):
                ef = em + es * float(ls) + eh * float(lh)
            else:
                x = float(ls) / 6.0e0
                ef = em + es * 6.0e0 + eh * float(lh)
                ef = ef + c * math.log(x)
            #endif
        #endif

        ip = ii + 1
        jp = rna.ibsp[ip]

        if ( jp != 0 ):
            e5 = EDANGLE(rna.iseq,ip,jp,ii,n)
        #endif

        if ( ii > 1 ):

            ip = ii - 1
            jp = rna.ibsp[ip]

            if ( jp != 0 ):
                ed = EDANGLE(rna.iseq,ip,jp,ii,n)
                e5 = min(e5,ed)
            #endif
        #endif

        ip = jj - 1
        jp = rna.ibsp[ip]

        if ( jp != 0 ):
            e3 = EDANGLE(rna.iseq,ip,jp,jj,n)
        #endif

        if ( jj < n ):

            ip = jj + 1
            jp = rna.ibsp[ip]

            if ( jp != 0 ):
                ed = EDANGLE(rna.iseq,ip,jp,jj,n)
                e3 = min(e3,ed)
            #endif
        #endif

        if ( mh == 2 and ms == 0 ):

            hs = rna.iseq[ii+1]
            js = rna.iseq[jj-1]

            ef += eaup[hs][js]
        
        #endif

        ef = ef + e3 + e4 + e5

    #endif

    dg = float(ef) - float(ei)

    return dg
