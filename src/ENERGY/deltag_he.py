"""
Module: DELTAG_HE (RNA,II,JJ)

Description: Computes the difference in free energy of an RNA loop due to
             extension of the helix below the ii-jj base pair using the
             empirical INN model.

Arguments:

          RNA - Class containing information about the RNA fold and
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

from rnavar import em,eh,es,eaup,beta

def DELTAG_HE(rna,ii,jj):

    # INTEGER
    # ii,jj
    # i,j,k,n,ip,jp,hs,js
    # nh,ns,mh,ms,iloop,indx
    # hs replaced is from FORTRAN code

    # REAL
    # e1,e3,e5,ei,ef
    # x,c,ed

    indx = rna.link[jj]

    i = rna.loop[indx]
    j = rna.ibsp[i]
    n = rna.n

    nh = rna.nhlx[indx]
    ns = rna.nsgl[indx]

    c = 1.75e0 / float(beta)

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
    ms = ns - 2

    #=== INITIAL ENERGY ===#

    if ( nh == 1 and iloop == 1 ):
        
        ei = EHAIR(rna.iseq,i,j,n)

    elif ( nh == 2 and iloop == 1 ):

        ip = i + 1
        while ( rna.ibsp[ip] == 0 ):
            ip += 1
        #endwhile

        jp = rna.ibsp[ip]

        ei = EBULGE(rna.iseq,i,j,ip,jp,n)

    else:

        ei = 0.0e0
        e5 = 0.0e0
        e3 = 0.0e0

        hs = rna.iseq[ii]
        js = rna.iseq[jj]

        if ( iloop == 1 ):

            if ( ns <= 6 ):
                ei = em + es * float(ns) + eh * REAL(nh)
            else:
                x = float(ns) / 6.0e0
                ei = em + es * 6.0e0 + eh * float(nh)
                ei = ei + c * math.log(x)
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
                #endif
            #endif

            if ( ii > 3 ) and ( rna.ibsp[ii-2] == 0 ):

                ip = ii - 3
                jp = rna.ibsp[ip]

                if ( jp != 0 ):
                    ed = EDANGLE(rna.iseq,ip,jp,ii-2,n)
                    ei += ed
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
                #endif
            #endif

            if ( jj < n-2 ) and ( rna.ibsp[jj+2] == 0 ):

                ip = jj + 3
                jp = rna.ibsp[ip]

                if ( jp != 0 ):
                    ed = EDANGLE(rna.iseq,ip,jp,jj+2,n)
                    ei += ed
                #endif
            #endif
        #endif

        ei = ei + e3 + e5 + eaup[hs][js]

    #endif


    #=== FINAL ENERGY ===#

    ef = ESTACK(rna.iseq,ii-1,jj+1,ii,jj,n)

    if ( mh == 1 and iloop == 1 ):

        e1 = EHAIR(rna.iseq,i+1,j-1,n)

    elif ( mh == 2 and iloop == 1 ):

        if ( ms == 0 ):

            e1 = ESTACK(rna.iseq,ii-2,jj+2,ii-1,jj+1,n)

        else:

            ip = i + 1
            while ( rna.ibsp[ip] == 0 ):
                ip += 1
            #endwhile

            jp = rna.ibsp[ip]

            if ( j == ii ):
                e1 = EBULGE(rna.iseq,i+1,j-1,ip,jp,n)
            else:
                e1 = EBULGE(rna.iseq,i,j,ip-1,jp+1,n)
            #endif
        #endif

    else:

        e1 = 0.0e0
        e5 = 0.0e0
        e3 = 0.0e0


        hs = rna.iseq[ii-1]
        js = rna.iseq[jj+1]

        if ( iloop == 1 ):

            if ( ms <= 6 ):
                e1 = em + es * float(ms) + eh * float(mh)
            else:
                x = float(ms) / 6.0e0
                e1 = em + es * 6.0e0 + eh * float(mh)
                e1 = e1 + c * math.log(x)
            #endif
        #endif

        if ( ii > 2 ) and ( rna.ibsp[ii-2] == 0 ):

            e5 = EDANGLE(rna.iseq,ii-1,jj+1,ii-2,n)

            if ( ii > 3 ):

                ip = ii - 3
                jp = rna.ibsp[ip]

                if ( jp != 0 ):
                    ed = EDANGLE(rna.iseq,ip,jp,ii-2,n)
                    e5 = min(e5,ed)
                #endif
            #endif
        #endif

        if ( jj < n-1 ) and ( rna.ibsp[jj+2] == 0 ):

            e3 = EDANGLE(rna.iseq,ii-1,jj+1,jj+2,n)

            if ( jj < n-2 ):

                ip = jj + 3
                jp = rna.ibsp[ip]

                if ( jp != 0 ):
                    ed = EDANGLE(rna.iseq,ip,jp,jj+2,n,ed)
                    e3 = min(e3,ed)
                #endif
            #endif
        #endif

        e1 = e1 + e3 + e5 + eaup[hs][js]

    #endif

    ef += e1

    dg = float(ef) - float(ei)

    return dg




