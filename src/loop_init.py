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

    #=== Find loops ===#

    nl = 1
    rna.loop[0] = n

    for i in range(0,n):

        j = rna.ibsp[i]

        if ( j > i ):

            ip = i + 1
            jp = j - 1

            if ( rna.ibsp[ip] != jp ):
                nl = nl + 1
                rna.loop(nl) = i

    rna.nl = nl

    #=== Make links ===#

    for i in range(0,nl):

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

        kp = ks

        while ( kp <= ke ):

            if ( rna.ibsp[kp] > kp ): nh += 1
            if ( rna.ibsp[kp] == 0 ): ns += 1

            if ( rna.ibsp[kp] > kp ):
                kp = rna.ibsp[kp]
                rna.link[kp] = i
            else:
                kp += 1

        rna.nhlx[i] = nh
        rna.nsgl[i] = ns

    #=== Compute size of partial sum table ===#

    nsum = 2
    while ( nsum < nl ):
        nsum *= 2
    rna.nsum = nsum

    #=== Compute reactions for loops ===#
    
    for i in range(0,nl):
        rna.LOOP_REAC(i)

    return rna