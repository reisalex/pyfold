"""
Module: ESTUC (ISEQ,IBSP,N)

Description: Computes the energy of an RNA secondary structure using the
             empirical INN model.

Method: Uses the MFOLD 3.0 energy function for RNA given by:

        E = SUM E_loop + E_stack

Arguments:

          ISEQ - Array of length N containing the sequence
                 in numerical code (A=1,C=2,G=3,U=4)
          IBSP - Array of dimension (N) containing the information
                 on base pairs in the RNA fold.
                 IBSP(i) = j  [i base pairs with j]
                 IBSP(i) = -1 [i is single stranded]
             N - Number of nucleotides in the sequence.

History:
Version     Date            Comment
--------    -------         --------------------
            09/28/2017      Original Code

Dependencies:

Author(s): Alex Reis
           Copyright (c) 2017 (Please refer to LICENCE)
"""

def ESTRUC(iseq,ibsp,n):

    # INTEGER
    # iseq(n),ibsp(n)
    # i,j,ip,jp
    # il,nl,loop(n)

    # REAL
    # e,el,es

    loop = [0 for _ in range(n)]

    e  = 0.0e0
    el = 0.0e0
    es = 0.0e0

    #=== Find Loops ===#

    nl = 1
    loop[0] = n - 1

    for i in xrange(0,n):

        j = ibsp[i]

        if ( j > i ):

            ip = i + 1
            jp = j - 1

            if ( ibsp[ip] != jp ):
                loop[nl] = i
                nl += 1

    #=== Compute Energy ===#

    for il in range(0,nl):

        i = loop[il]
        j = ibsp[i]

        if ( i == n - 1 ): j = 0

        #=== Loop Energy ===#

        el = ELOOP(iseq,ibsp,i,j,n)

        e += el

        #=== Stacking Energy ===#

        ip = i - 1
        jp = j + 1

        if ( ip < 0   ): continue
        if ( jp > n-1 ): continue
        if ( i  > j   ): continue

        while ( ibsp[ip] == jp ):

            es = ESTACK(iseq,ip,jp,i,j,n)

            e += es

            i = ip
            j = jp

            ip -= 1
            jp += 1

            if ( ip < 0   ): break
            if ( jp > n-1 ): break

    return e

