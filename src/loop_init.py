"""
Subroutine: LOOP_INIT (RNA)

Description: Initializes the data structures (loop elements) required
             for RNA kinetics.

Arguments:
        
        RNA - Class structure containing information on the
                RNA secondary structure and possible reactions

History:
Version     Date            Comment
--------    -------         --------------------
            09/28/2017      Original Code

Dependencies:

Author(s): Alex Reis
             Copyright (c) 2017 (Please refer to LICENCE)
"""

def LOOP_INIT(rna):

    # VARIABLES
    
    # INTEGERS
    # i,j,n,nl,ns,nh
    # ip,jp,kp,ks,ke
    # nsum

    #=== Initialize RNA Data ===#

    n = rna.n
    rna.CLEAR_LOOPS()

    #=== Find loops ===#

    nl = 1
    rna.loop[1] = n

    for i in range(1,n):

        j = rna.ibsp[i]

        if ( j > i ):

            ip = i + 1
            jp = j - 1

            if ( rna.ibsp[ip] != jp ):
                nl += 1
                rna.loop[nl] = i

    rna.nl = nl

    #=== Make links ===#

    for i in range(1,nl):

        ip = rna.loop[i]
        jp = rna.ibsp[ip]

        if ( ip == n ): jp = 1
        
        rna.link[ip] = i

        if ( ip < jp ):
            ks = ip + 1
            ke = jp
            nh = 1
            ns = 0
        else:
            ks = jp
            ke = ip
            nh = 0
            ns = 0

        ins = 0
        kp = ks

        # modified to collect loop asymmetry
        while ( kp <= ke ):

            # unpaired nt
            if ( rna.ibsp[kp] == -1 ):
                ins += 1
                kp += 1

            # new helix in loop
            elif ( ibsp[kp] > kp ):
                rna.lns[i][nh] = ins
                ns += ins
                ins = 0
                nh += 1

                # track helix ip-jp in MBLs
                rna.htrack[i][kp] = nh

                # skip to closing bp of current nt
                kp = rna.ibsp[kp]
                rna.link[kp] = i

        rna.lns[i][0]  = ins # for-loop simplicity
        rna.lns[i][nh] = ins
        ns += ins

        rna.nhlx[i] = nh
        rna.nsgl[i] = ns

    # Compute size of partial sum table

    nsum = 2
    while ( nsum < nl ):
        nsum *= 2
    rna.nsum = nsum

    # Compute reactions for loops
    
    for i in range(1,nl):
        rna.LOOP_REAC(i)

    return rna