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
                 IBSP(i) = j  [i base pairs with j]
                 IBSP(i) = -1 [i is single stranded]
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
from rnavar import em,eh,es,beta

def ELOOP(iseq,ibsp,i,j,n):

    # INTEGER
    # i,j,n,nh,ns,tns,indx
    # iseq(n),ibsp(n),nsl(n)

    # REAL
    # x,c,ed,e3,e5,el,asym

    el = 0.0e0

    nh = 0
    ns = 0
    ins = 0
    lns = []

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

        # unpaired nt
        if ( ibsp[k] == -1 ):
            ins += 1

        # new helix in loop
        elif ( ibsp[k] >  k ):
            lns[nh] = ins
            ns += ins
            ins = 0
            nh += 1

            # skip to closing bp of current nt
            if ( iloop == 0 or k != ks ):
                k = ibsp[k]

        k += 1

    lns[0]  = ins # for-loop simplicity
    lns[nh] = ins
    ns += ins

    #=== Compute Loop Energy ===#

    if ( nh == 1 and iloop == 1 ):

        el = EHAIR(iseq,i,j,n)

    elif ( nh == 2 and iloop == 1 ):

        ip = i + 1
        while ( ibsp[ip] == -1 ):
            ip += 1

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

            if ( k >= 0 and k <= n-1 ) and ( ibsp[k] == -1 ):
                e3 = EDANGLE(iseq,ip,jp,k,n)

            if ( params.MBLmodel == 2 ):

                # calculate average asymmetry of MBL
                for indx in xrange(1,nh+1):
                    asym += float( abs( lns[indx] - lns[indx-1] ) )
                asym /= float(nh)
                asym = min(2.0,asym)

                el  = params.MBLinit[0]               # a
                el += params.MBLinit[1] * float(asym) # b
                el += params.MBLinit[2] * float(nh)   # c

                if (nh == 3) and (ns < 2):
                    el += params.MBLinit[4] # dG_strain

            elif ( ns < 6 ):
                el  = params.MBLinit[0]               # a, (em)
                el += params.MBLinit[1] * float(ns)   # b, (es)
                el += params.MBLinit[2] * float(nh)   # c, (eh)

            else:
                el  = params.MBLinit[0]
                el += params.MBLinit[1] * 6.0e0
                x   = float(ns) / 6.0e0
                el += c * math.log(x)
                el += params.MBLinit[2] * float(nh)

        while ( k <= ke ):

            if ( ibsp[k] != -1 ):

                ip = k
                jp = ibsp[k]

                kp = ip - 1

                e5 = 0.0e0

                if ( kp >= 0 and kp <= n-1 ) and ( ibsp[kp] == -1 ):
                    e5 = EDANGLE(iseq,ip,jp,kp,n)

                # pick 3' (H1) or 5' (H2) dangle if same nt
                if ( ilast == kp ):
                    ed = min(e3,e5)
                else:
                    ed = e3 + e5

                el = el + ed + params.dG_AUP[iseq[ip]][iseq[jp]]

                kp = jp + 1
                ilast = kp

                e3 = 0.0e0

                if ( kp >= 0 and kp <= n-1 ) and ( ibsp[kp] == -1 ):
                    e3 = EDANGLE(iseq,ip,jp,kp,n)

                if ( k != ke ):
                    k = ibsp[k]
            k += 1

        if ( iloop == 0 ):
            el += e3

    return el