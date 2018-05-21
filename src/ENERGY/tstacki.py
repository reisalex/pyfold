'''
Subroutine: TSTACKI (LIST,ES)

Purpose: Performs a table lookup for the stacking interaction of
         the two nucleotides positioned over a closing basepair in
         an internal loop.

Method: Uses the MFOLD 3.0 energy function table for RNA @ T=37.

Arguments:

         LIST - Array of length 4 containing the nucleotides in
                numerical code (A=1,C=2,G=3,U=4) for the
                following locations:

                5' (1) A X (3) 3' INTERNAL LOOP
                3' (2) U Y (4) 5' INTERNAL LOOP

                where LIST(1) = letter code for position 1 etc.

           ES - (OUTPUT) MFOLD 3.0 stacking energy of the sequence
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

def TSTACKI(ilist,es):

    # INTEGERS
    # ilist(4),i,i1,i2,i3,i4

    # REAL
    # es

    i1 = ilist[0]
    i2 = ilist[1]
    i3 = ilist[2]
    i4 = ilist[3]

    es += params.dG_stacki[i1,i2,i3,i4]

    return es