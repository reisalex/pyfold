"""
Module: DELTAG_HI (RNA,I,J)

Description: Computes the difference in free energy of an RNA helix due to a
             deletion of the internal i-j base pair using the INN model.

Arguments:

            RNA - Class containing information about the RNA fold and
                 the RNA sequence.
             I  - Nucleotide position of the 5' nucleotide.
             J  - Nucleotide position of the 3' nucleotide.

History:
Version     Date            Comment
--------    -------         --------------------
            09/28/2017      Original Code

Dependencies:

Author(s): Alex Reis
           Copyright (c) 2017 (Please refer to LICENCE)
"""

def DELTAG_HI(rna,i,j):

    # INTEGER
    # i,j,n

    # REAL
    # ei,ef,es

    # DOUBLE PRECISION
    # dg

    n = rna.n

    #=== Initial Energy ===#

    ei = ESTACK(rna.iseq,i-1,j+1,i,j,n)
    es = ESTACK(rna.iseq,i,j,i+1,j-1,n)

    ei += es

    #=== Final Energy ===#

    ef = EBULGE(rna.iseq,i-1,j+1,i+1,i-1,n)

    dg = float(ef) - float(ei)

    return dg