'''
Subroutine: TINT11 (LIST,EB)

Purpose: Performs a table lookup for the internal energy for a mismatch
         pair between two basepairs in a helix.

Method: Uses the MFOLD 3.0 energy function table for RNA @ T=37.

Arguments:

         LIST - Array of length 4 containing the nucleotides in
                numerical code (A=1,C=2,G=3,U=4) for the
                following locations:

                        (3)         
                5' (1) A . X (5) 3' 
                3' (2) U . Y (6) 5' 
                        (4)         

                where LIST(1) = letter code for position 1 etc.

           EB - (OUTPUT) MFOLD 3.0 internal loop energy of the sequence
                provided in LIST.

History:

Version    Date         Comment
--------   ----------   -----------------------
           01/01/2015   Original Code

Dependencies:

Modules -
Functions -
Subroutines -

Author(s): Alex Reis
           Copyright (c) 2015 (Please Refer to LICENCE)
'''

def TINT11(ilist,eb):

    # INTEGERS
    # ilist(6),i,j,i1,i2,i3,i4,i5,i6

    # REAL
    # eb
    
    i1 = ilist[0]
    i2 = ilist[1]
    i3 = ilist[2]
    i4 = ilist[3]
    i5 = ilist[4]
    i6 = ilist[5]

    eb += params.dG_int11[i1,i2,i3,i4,i5,i6]

    return eb