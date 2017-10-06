"""
Function: ELOOP (ISEQ,IBSP,I,J,N)

Description: Computes the energy of an RNA loop using the empirical INN model.

Method: Uses the MFOLD 3.0 energy function for RNA given by:
           
           EL = EHAIR                --> FOR k = 1
           EL = EBULGE               --> FOR k = 2
                                     --> FOR k > 2

           EL = a + b*NS + c*NH + GS               IF NS <= 6
              = a + b *6 + c*NH + GS + d*LOG(NS/6) IF NS  > 6

           where k = number of helices in loop
                 a = Multi-loop Penality
                 b = Single Stranded Penalty
                 c = Helix Penalty
                 d = 1.75 KT
                 GS= Stacking Energy of Dangling Bases

Arguments:

          ISEQ - Array of length N containing the sequence
                 in numerical code (A=1,C=2,G=3,U=4)
          IBSP - Array of dimension (N) containing the information
                 on base pairs in the RNA fold.
                 IBSP(i) = j [i base pairs with j]
                 IBSP(i) = 0 [i is single stranded]
             I - Nucleotide position of the 5' most nucleotide.
             J - Nucleotide position of the 3' most nucleotide.
             N - Number of nucleotides in the sequence.

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

def ELOOP(iseq,ibsp,i,j,n):

    # INTEGER
    # i,j,n
    # iseq(n),ibsp(n)

    # REAL
    # x,c,ed,e3,e5,el

    el = 0.0e0

    nh = 0
    ns = 0

    c = 1.750e0 / float(beta)

    #=== Internal Loop iloop = 1 ===#
    #=== External Loop iloop = 0 ===#

    if ( i < j ): iloop = 1
    if ( i > j ): iloop = 0

    ks = min(i,j)
    ke = max(i,j)

    #=== Count number of helices in loop ===#

    k = ks

    while ( k <= ke ):

        if ( ibsp[k] == 0 ): ns += 1
        if ( ibsp[k] >  k ): nh += 1

        if ( ibsp[k] > k ) and ( iloop == 0 or k != ks ):
            k = ibsp[k]
        #endif

        k += 1

    #endwhile

    #=== Compute Loop Energy ===#

    if ( nh == 1 and iloop == 1 ):

        el = EHAIR(iseq,i,j,n)

    elif ( nh == 2 and iloop == 1 ):

        ip = i + 1
        while ( ibsp[ip] == 0 ):
            ip += 1
        #endwhile

        jp = ibsp[ip]

        el = EBULGE(iseq,i,j,ip,jp,n)

    else:

        e3 = 0.0e0

        if ( iloop == 0 ):

            k = ks
            ilast = 0

        elif ( iloop == 1 ):

            ip = i
            jp = ibsp[i]

            k = ip + 1
            ilast = k

            if ( k >= 1 and k <= n ) and ( ibsp[k] == 0 ):
                e3 = EDANGLE(iseq,ip,jp,k,n)
            #endif

            if ( ns <= 6 ):
                el = em + es * float(ns) + eh * float(nh)
            else:
                x = float(ns) / 6.0e0
                el = em + es * 6.0e0 + eh * float(nh)
                el = el + c * math.log(x)
            #endif

            while ( k <= ke ):

                if ( ibsp[k] != 0 ):

                    ip = k
                    jp = ibsp[k]

                    kp = ip - 1

                    e5 = 0.0e0

                    if ( kp >= 1 and kp <= n ) and ( ibsp[kp] == 0 ):
                        e5 = EDANGLE(iseq,ip,jp,kp,n)
                    #endif

                    if ( ilast == kp ):
                        ed = min(e3,e5)
                    else:
                        ed = e3 + e5
                    #endif

                    el = el + ed + eaup[iseq[ip]][iseq[jp]]

                    kp = jp + 1
                    ilast = kp

                    e3 = 0.0e0

                    if ( kp >= 1 and kp <= n ) and ( ibsp[kp] == 0 ):
                        e3 = EDANGLE(iseq,ip,jp,kp,n)
                    #endif

                    if ( k != ke ):
                        k = ibsp[k]
                    #endif
                #endif

                k += 1

            #endwhile

            if ( iloop == 0 ):
                el += e3
            #endif
        #endif

        return el