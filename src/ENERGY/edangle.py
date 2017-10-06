"""
Function: EDANGLE (ISEQ,I,J,K,N)

Description: Computes the energy of a dangling nucleotide over a basepair
             using the empirical INN model.

Method: Uses the MFOLD 3.0 energy function for RNA @ T=37 given by:

        ED = E_dangle

              5' (I) X       3'
              3' (J) Y Z (K) 5'


Arguments:

         ISEQ - Array of length N containing the sequence
                in numerical code (A=1,C=2,G=3,U=4)
            I - Nucleotide position of the basepair 5'.
            J - Nucleotide position of the basepair 3'.
            K - Nucleotide position of the dangling nucleotide.
            N - Number of nucleotides in the sequence.

History:
Version     Date            Comment
--------    -------         --------------------
            09/28/2017      Original Code

Dependencies:

Author(s): Alex Reis
           Copyright (c) 2017 (Please refer to LICENCE)
"""

def EDANGLE(iseq,i,j,k,n):

    # INTEGER
    # i,j,k,n
    # iseq(n)
    # list(3); now called ilist(3)

    # REAL
    # ed

    ed = 0.0e0

    if ( k < 1 ) return ed
    if ( k > n ) return ed

    if ( k == i+1 ):

        # 5' (i) A X (k) 3'
        # 3' (j) U       5'

        ilist[0] = iseq[i]
        ilist[1] = iseq[j]
        ilist[2] = iseq[k]

        ed = TDANGLE3(ilist)

    #endif

    if ( k == j-1 ):

        # 5' (i) A       3'
        # 3' (j) U X (k) 5'        

        ilist[0] = iseq[i]
        ilist[1] = iseq[j]
        ilist[2] = iseq[k]

        ed = TDANGLE5(ilist)

    #endif

    if ( k == j+1 ):

        # 5'       A (i) 3'
        # 3' (k) X U (j) 5'

        ilist[0] = iseq[j]
        ilist[1] = iseq[i]
        ilist[2] = iseq[k]

        ed = TDANGLE3(ilist)

    #endif

    if ( k == i-1 ):

        # 5' (k) X A (i) 3'
        # 3'       U (j) 5'

        ilist[0] = iseq[j]
        ilist[1] = iseq[i]
        ilist[2] = iseq[k]

        ed = TDANGLE5(ilist)

    #endif

    return ed