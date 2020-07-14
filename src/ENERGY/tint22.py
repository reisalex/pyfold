'''
Subroutine: TINT22 (LIST,EB)

Purpose: Performs a table lookup for the internal energy for two mismatch
         pairs between two basepairs in a helix.

Method: Uses the MFOLD 3.0 energy function table for RNA @ T=37.

Arguments:

         LIST - Array of length 4 containing the nucleotides in
                numerical code (A=1,C=2,G=3,U=4) for the
                following locations:

                        (3)(5)        
                5' (1) A .  . X (7) 3'
                3' (2) U .  . Y (8) 5'
                        (4)(6)             

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

Author(s): Eric Dykeman
           Copyright (c) 2015 (Please Refer to LICENCE)
'''

def TINT22(ilist,eb):

    # INTEGERS
    # ilist(8),i1,i2,i3,i4,i5,i6,i7,i8

    # REAL
    # eb

    i1 = ilist[0]
    i2 = ilist[1]
    i3 = ilist[2]
    i4 = ilist[3]
    i5 = ilist[4]
    i6 = ilist[5]
    i7 = ilist[6]
    i8 = ilist[7]

    eb += params.dG_int22[i1,i2,i3,i4,i5,i6,i7,i8]

    return eb