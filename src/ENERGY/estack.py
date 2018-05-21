"""
Function: ESTACK (ISEQ,I,J,IP,JP,N)

Description: Computes the energy of a helix stacking between two bp using
             the empirical INN energy model.

Method: Uses the MFOLD 3.0 energy function for RNA given by:

        ES = E_stack

              5' (I) X W (IP) 3'
              3' (J) Y Z (JP) 5'

              NOTE: I < IP and JP < J

Arguments:

         ISEQ - Array of length N containing the sequence
                in numerical code (A=1,C=2,G=3,U=4)
            I - Nucleotide position of the first basepair 5'.
            J - Nucleotide position of the first basepair 3'.
           IP - Nucleotide position of the second basepair 5'.
           JP - Nucleotide position of the second basepair 3'.
            N - Number of nucleotides in the sequence.

History:
Version     Date            Comment
--------    -------         --------------------
            09/28/2017      Original Code

Dependencies:

Author(s): Alex Reis
           Copyright (c) 2017 (Please refer to LICENCE)
"""

def ESTACK(iseq,i,j,ip,jp,n):

    # INTEGER
    # i,j,ip,jp,n
    # iseq(n)
    # ilist(4)

    # REAL
    # es

    es = 0.0e0

    # 5' (i) A X (ip) 3'
    # 3' (j) U Y (jp) 5'

    ilist[0] = iseq[i]
    ilist[1] = iseq[j]
    ilist[2] = iseq[ip]
    ilist[3] = iseq[jp]

    es = TSTACK(ilist,es)

    return es
