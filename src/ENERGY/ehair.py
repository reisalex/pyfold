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
from rnavar import eau,beta

def EHAIR(iseq,i,j,n):

    # INTEGER
    # i,j,n
    # iseq(n)
    # ilist(6)
    # k,ic,nl

    # REAL
    # x,c,elh(30)
    # eh

    elh = [0.00e0,0.00e0,5.70e0,5.60e0,5.60e0,
           5.40e0,5.90e0,5.60e0,6.40e0,6.50e0,
           6.60e0,6.70e0,6.80e0,6.90e0,6.90e0,
           7.00e0,7.10e0,7.10e0,7.20e0,7.20e0,
           7.30e0,7.30e0,7.40e0,7.40e0,7.50e0,
           7.50e0,7.50e0,7.60e0,7.60e0,7.70e0]

    eh = 0.0e0

    nl = j - i - 1

    c = 1.750e0 / float(beta)

    #=== TERM 1 --> Entropic Term ===#

    if ( nl > 30 ):

        x = float(nl) / 30.0e0
        x = c * math.log(x)

        eh = elh[30] + x

    else:

        eh = elh[nl]

    #endif

    #=== TERM 2 --> Stacking Energy ===#

    if ( nl > 3 ):

        # 5' (i) A X (i+1) LOOP
        # 3' (j) U Y (j+1) LOOP

        ilist[0] = iseq[i]
        ilist[1] = iseq[j]
        ilist[2] = iseq[i+1]
        ilist[3] = iseq[j-1]

        eh = TSTACKH(ilist,eh)

    #endif

    #=== TERM 3 ---> Bonuses ===#

    #=== Tetra-loop Bonus ===#

    if ( nl == 4 ):

        ilist[0] = iseq[i]
        ilist[1] = iseq[i+1]
        ilist[2] = iseq[i+2]
        ilist[3] = iseq[i+3]
        ilist[4] = iseq[i+4]
        ilist[5] = iseq[i+5]

        eh = TLOOP(ilist,eh)

    #endif

    #=== GGG Hairpin Bonus ===#

    if ( iseq[i] == 3 and iseq[j] == 4 ):

        ic = 0

        for k in range(max(1,i-2),i):
            if ( iseq[k] == 3 ): ic += 1
        #endfor

        if ( ic == 3 ): eh -= 2.20e0

    #endif

    #=== TERM 4 --> Penalties ===#

    #=== Poly C Penalty ===#

    ic = 0

    for k in range(i+1,j-1):
        if ( iseq[k] == 2 ): ic += 1
    #endfor

    if ( ic == nl ):
        if ( nl == 3 ):
            eh += 1.40e0
        else:
            eh += 1.60e0
            eh += 0.30e0 * float(ic)
        #endif
    #endif

    #=== A-U / G-U closing a Tri-loop ===#

    if ( nl == 3 ):

        if ( iseq[i] == 4 ): eh += eau
        if ( iseq[j] == 4 ): eh += eau

    #endif

    return eh
        
            


