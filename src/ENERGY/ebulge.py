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
    # x,c,f(4)

    eb = 0.0e0

    n1 = ip - i - 1 # number of nt on 5' side of bulge/int loop
    n2 = j - jp - 1 # number of nt on 3' side of bulge/int loop

    nt = n1 + n2 # loop size
    na = abs(n1-n2) # asymmetry

    imin = min(n1,n2)
    imax = max(n1,n2)

    c = 1.750e0 / float(beta)

    #=== Get Bulge Type ===#

    ibul = 6
    if   ( n1 == 0 and n2 == 1 ): ibul = 0
    elif ( n1 == 1 and n2 == 0 ): ibul = 0
    elif ( n1 == 0 and n2 >= 2 ): ibul = 1
    elif ( n1 >= 2 and n2 == 0 ): ibul = 1
    elif ( n1 == 1 and n2 == 1 ): ibul = 2
    elif ( n1 == 1 and n2 == 2 ): ibul = 3
    elif ( n1 == 2 and n2 == 1 ): ibul = 4
    elif ( n1 == 2 and n2 == 2 ): ibul = 5

    #=== TERM 1 --> Entropic Term ===#

    if ( nt > 30 ):

        x = float(nt) / 30.0e0
        x = c * math.log(x)

        if ( ibul <= 1 ): eb = params.dG_bulge[30] + x
        if ( ibul == 6 ): eb = params.dG_iloop[30] + x

    else:

        if ( ibul <= 1 ): eb = params.dG_bulge[nt]
        if ( ibul == 6 ): eb = params.dG_iloop[nt]

    #endif

    #=== TERM 2 & 3 --> Stacking Energy ===#

    if ibul == 0:

        # 5' (i) A . X (ip) 3'
        # 3' (j) U   Y (jp) 5'

        ilist[0] = iseq[i]
        ilist[1] = iseq[j]
        ilist[2] = iseq[ip]
        ilist[3] = iseq[jp]

        eb = TSTACK(ilist,eb)

    elif ibul == 1:

        # 5' (i) A .. X (ip) 3'
        # 3' (j) U    Y (jp) 5'

        #=== Closing A-U / G-U Penalty ===#

        if ( iseq[i]  == 4 ): eb += params.dG_AU
        if ( iseq[j]  == 4 ): eb += params.dG_AU
        if ( iseq[ip] == 4 ): eb += params.dG_AU
        if ( iseq[jp] == 4 ): eb += params.dG_AU

    elif ibul == 2:

        # 5' (i) A . X (ip) 3'
        # 3' (j) U . Y (jp) 5'

        ilist[0] = iseq[i]
        ilist[1] = iseq[j]
        ilist[2] = iseq[i+1]
        ilist[3] = iseq[j-1]
        ilist[4] = iseq[ip]
        ilist[5] = iseq[jp]

        eb = TINT11(ilist,eb)

    elif ibul == 3:

        # 5' (i) A .   X (ip) 3'
        # 3' (j) U . . Y (jp) 5'

        ilist[0] = iseq[i]
        ilist[1] = iseq[j]
        ilist[2] = iseq[i+1]
        ilist[3] = iseq[j-1]
        ilist[4] = iseq[j-2]
        ilist[5] = iseq[ip]
        ilist[6] = iseq[jp]

        eb = TINT12(ilist,eb)

    elif ibul == 4:

        # 5' (i) A . . X (ip) 3'
        # 3' (j) U   . Y (jp) 5'

        ilist[0] = iseq[jp]
        ilist[1] = iseq[ip]
        ilist[2] = iseq[jp+1]
        ilist[3] = iseq[ip-1]
        ilist[4] = iseq[ip-2]
        ilist[5] = iseq[j]
        ilist[6] = iseq[i]

        eb = TINT12(ilist,eb)

    elif ibul == 5:

        # 5' (i) A . . X (ip) 3'
        # 3' (j) U . . Y (jp) 5'

        ilist[0] = iseq[i]
        ilist[1] = iseq[j]
        ilist[2] = iseq[i+1]
        ilist[3] = iseq[j-1]
        ilist[4] = iseq[ip-1]
        ilist[5] = iseq[jp+1]
        ilist[6] = iseq[ip]
        ilist[7] = iseq[jp]

        eb = TINT22(ilist,eb)

    elif ibul == 6:

        # 5' (i) A X .. G (ip) 3'
        # 3' (j) U Y .. C (jp) 5'

        ilist[0] = iseq[i]
        ilist[1] = iseq[j]
        ilist[2] = iseq[i+1]
        ilist[3] = iseq[j-1]

        #=== GAIL Rule ===#

        if ( imin == 1 and imax > 2 ):

            ilist[2] = 0
            ilist[3] = 0

        eb = TSTACKI(ilist,eb)

        # 5' (i) A .. X G (ip) 3'
        # 3' (j) U .. Y C (jp) 5'

        ilist[0] = iseq[jp]
        ilist[1] = iseq[ip]
        ilist[2] = iseq[jp+1]
        ilist[3] = iseq[ip-1]

        #=== GAIL Rule ===#

        if ( imin == 1 and imax > 2 ):

            ilist[2] = 0
            ilist[3] = 0

        eb = TSTACKI(ilist,eb)

    else:

        raise Exception('ERROR: ibul out of range [0-6]')

    #=== PART 4 ---> Asymmetry Penalty ===#

    if ( ibul == 6 ):

        k = min(5,n1,n2)

        x = float(na) * params.dG_asym[k]

        x = min(x,params.dG_maxasym)

        eb += x

    return eb