"""
Function: EBULGE (ISEQ,I,J,IP,JP,N)

Description: Computes the energy of an RNA bulge with two helices using the
             empirical INN energy model.

Method: Uses the MFOLD 3.0 energy function for RNA given by:

        EB = E_entropic + E_stack + E_stack + E_asymmetry

              5' (I) X ... W (IP) 3'
              3' (J) Y ... Z (JP) 5'

              NOTE: I < IP and JP < J

Arguments:

         ISEQ - Array of length N containing the sequence
                in numerical code (A=1,C=2,G=3,U=4)
            I - Nucleotide position of the starting basepair 5'.
            J - Nucleotide position of the starting basepair 3'.
           IP - Nucleotide position of the ending basepair 3'.
           JP - Nucleotide position of the ending basepair 5'.
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
from rnavar import eau,beta

def EBULGE(iseq,i,j,ip,jp,n):

    # INTEGER
    # i,j,ip,jp,n
    # iseq(n)
    # k,ibul,imin,imax
    # n1,n2,nt,na,ilist(8)

    # REAL
    # eb
    # x,c,f(4),elb(30),eli(30)
    

    #=== MOVE THESE DATA TABLES TO DATABASE FOR ENERGY PARAMS ===#

    f = [0.50e0,0.50e0,0.50e0,0.50e0]

    elb = [3.80e0,2.80e0,3.20e0,3.60e0,4.00e0,
           4.40e0,4.60e0,4.70e0,4.80e0,4.90e0,
           5.00e0,5.10e0,5.20e0,5.30e0,5.40e0,
           5.40e0,5.50e0,5.50e0,5.60e0,5.70e0,
           5.70e0,5.80e0,5.80e0,5.80e0,5.90e0,
           5.90e0,6.00e0,6.00e0,6.00e0,6.10e0]

    eli = [0.00e0,0.00e0,0.00e0,1.70e0,1.80e0,
           2.00e0,2.20e0,2.30e0,2.40e0,2.50e0,
           2.60e0,2.70e0,2.80e0,2.90e0,3.00e0,
           3.00e0,3.10e0,3.10e0,3.20e0,3.30e0,
           3.30e0,3.40e0,3.40e0,3.40e0,3.50e0,
           3.50e0,3.60e0,3.60e0,3.60e0,3.70e0]

    eb = 0.0e0

    n1 = ip - i - 1
    n2 = j - jp - 1

    nt = n1 + n2
    na = abs(n1-n2)

    imin = min(n1,n2)
    imax = max(n1,n2)

    x = 1.750e0 / float(beta)

    #=== Get Bulge Type ===#

    ibul = 6
    if ( n1 == 0 and n2 == 1 ): ibul = 0
    if ( n1 == 1 and n2 == 0 ): ibul = 0
    if ( n1 == 0 and n2 >= 2 ): ibul = 1
    if ( n1 >= 2 and n2 == 0 ): ibul = 1
    if ( n1 == 1 and n2 == 1 ): ibul = 2
    if ( n1 == 1 and n2 == 2 ): ibul = 3
    if ( n1 == 2 and n2 == 1 ): ibul = 4
    if ( n1 == 2 and n2 == 2 ): ibul = 5

    #=== TERM 1 --> Entropic Term ===#

    if ( nt > 30 ):

        x = float(nt) / 30.0e0
        x = c * math.log(x)

        if ( ibul <= 1 ): eb = elb[30] + x
        if ( ibul == 6 ): eb = eli[30] + x

    else:

        if ( ibul <= 1 ): eb = elb[nt]
        if ( ibul == 6 ): eb = eli[nt]

    #endif

    #=== TERM 2 & 3 --> Stacking Energy ===#

    if ibul == 0:

        # 5' (i) A . X (ip) 3'
        # 3' (j) U   Y (jp) 5'

        ilist[0] = iseq[i]
        ilist[1] = iseq[j]
        ilist[2] = iseq[ip]
        ilist[3] = iseq[jp]

        eb = TSTACK(ilist)

    elif ibul == 1:

        # 5' (i) A .. X (ip) 3'
        # 3' (j) U    Y (jp) 5'

        #=== Closing A-U / G-U Penalty ===#

        if ( iseq[i]  == 4 ): eb += eau
        if ( iseq[j]  == 4 ): eb += eau
        if ( iseq[ip] == 4 ): eb += eau
        if ( iseq[jp] == 4 ): eb += eau

    elif ibul == 2:

        # 5' (i) A . X (ip) 3'
        # 3' (j) U . Y (jp) 5'

        ilist[0] = iseq[i]
        ilist[1] = iseq[j]
        ilist[2] = iseq[i+1]
        ilist[3] = iseq[j-1]
        ilist[4] = iseq[ip]
        ilist[5] = iseq[jp]

        eb = TINT11(ilist)

    elif ibul == 3:

        # 5' (i) A .   X (ip) 3'
        # 3' (j) U . . Y (jp) 5'

