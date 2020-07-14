"""
Subroutine: EHAIR (ISEQ,I,J,N,EH)

Purpose: Computes the energy of an RNA hairpin turn using the
         empirical MFOLD 3.0 energy function.

Method: Uses the MFOLD 3.0 energy function for RNA @ T=37 given by:

        EH = E_entropic + E_stack + E_bonus + E_penalty

              5' (I) X ... loop 3'
              3' (J) Y ... loop 5'

              NOTE: I < J

Arguments:

         ISEQ - Array of length N containing the sequence
                in numerical code (A=1,C=2,G=3,U=4)
            I - Nucleotide position of the loop basepair 5'.
            J - Nucleotide position of the loop basepair 3'.
            N - Number of nucleotides in the sequence.
           EH - (OUTPUT) MFOLD 3.0 energy of the hairpin sequence.

History:
Version     Date            Comment
--------    -------         --------------------
            09/28/2017      Original Code

Dependencies:

Author(s): Alex Reis
           Copyright (c) 2017 (Please refer to LICENCE)
"""

import math
from rnavar import eau,beta

def EHAIR(iseq,i,j,n):

    # INTEGER
    # i,j,n
    # iseq(n)
    # ilist(8)
    # k,ic,nl

    # REAL
    # x,c
    # eh

    eh = 0.0e0

    nl = j - i - 1

    c = 1.750e0 / float(beta)

    #=== Tables of specific hairpin loops ===#

    #=== Triloop Bonus ===#

    if ( nl == 3 ):

        ilist[0] = iseq[i]
        ilist[1] = iseq[i+1]
        ilist[2] = iseq[i+2]
        ilist[3] = iseq[i+3]
        ilist[4] = iseq[i+4]

        eh = TLOOP(ilist,eh,nl)

        if eh != 0.0e0:
            return eh

    #=== Tetraloop Bonus ===#

    elif ( nl == 4 ):

        ilist[0] = iseq[i]
        ilist[1] = iseq[i+1]
        ilist[2] = iseq[i+2]
        ilist[3] = iseq[i+3]
        ilist[4] = iseq[i+4]
        ilist[5] = iseq[i+5]

        eh = TLOOP(ilist,eh,nl)

        if eh != 0.0e0:
            return eh

    #=== Hexaloop Bonus ===#

    elif ( nl == 6 ):

        ilist[0] = iseq[i]
        ilist[1] = iseq[i+1]
        ilist[2] = iseq[i+2]
        ilist[3] = iseq[i+3]
        ilist[4] = iseq[i+4]
        ilist[5] = iseq[i+5]
        ilist[6] = iseq[i+6]
        ilist[7] = iseq[i+7]

        eh = TLOOP(ilist,eh,nl)

        if eh != 0.0e0:
            return eh

    # If loop is not present in a table:

    #=== TERM 1 --> Entropic Term ===#

    if ( nl > 30 ):

        x = float(nl) / 30.0e0
        x = c * math.log(x)

        eh = params.dG_hloop[30] + x

    else:

        eh = params.dG_hloop[nl]

    #=== TERM 2 --> Stacking Energy ===#

    if ( nl > 3 ):

        # 5' (i) A X (i+1) LOOP
        # 3' (j) U Y (j+1) LOOP

        ilist[0] = iseq[i]
        ilist[1] = iseq[j]
        ilist[2] = iseq[i+1]
        ilist[3] = iseq[j-1]

        eh = TSTACKH(ilist,eh)

    #=== TERM 3 ---> Bonuses ===#

    #=== GGG Hairpin Bonus ===#

    if ( iseq[i] == 2 and iseq[j] == 3 ):

        ic = 0

        for k in xrange(max(0,i-2),i):
            if ( iseq[k] == 2 ): ic += 1

        if ( ic == 2 ): eh -= params.dG_bonuses[2]

    #=== TERM 4 --> Penalties ===#

    #=== Poly C Penalty ===#

    ic = 0

    for k in xrange(i+1,j):
        if ( iseq[k] == 1 ): ic += 1

    if ( ic == nl ):
        if ( nl == 3 ):
            eh += params.dG_bonuses[3]
        else:
            eh += params.dG_bonuses[5]
            eh += params.dG_bonuses[4] * float(ic)

    #=== A-U / G-U closing a Tri-loop ===#

    if ( nl == 3 ):

        if ( iseq[i] == 3 ): eh += eau
        if ( iseq[j] == 3 ): eh += eau

    # ^I think the AU penalty applies regardless of whether the
    # hairpin loop is a triloop or not

    return eh