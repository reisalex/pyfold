'''
Subroutine: TLOOP (LIST,EL)

Purpose: Performs a table lookup for the special bonus energies
         for various tetra-loops.

Method: Uses the MFOLD 3.0 energy function table for RNA @ T=37.

Arguments:

         LIST - Array of length N containing the nucleotides in
                numerical code (A=0,C=1,G=2,U=3) for the
                following locations:

                    1 2 3 4 5 6
                5'  A W X Y Z U   3'
                      L O O P

                where LIST(1) = letter code for position 1 etc.

           EL - (OUTPUT) MFOLD 3.0 tetra-loop bonus energy of the
                sequence provided in LIST.

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

def TLOOP(ilist,el,nl):

    # ARGUMENTS/VARIABLES
    # INTEGERS
    # ilist(8), nl, i
    # STRINGS
    # cl(8)
    # REAL INOUT
    # el

    # cwrk?


    # Triloops
    #      0 1 2 3 4      #
    #  5'  A W X Y U   3' #

    if nl == 3:

        for i in xrange(5):
            if   ilist[i] == 0: cl[i] = 'A'
            elif ilist[i] == 1: cl[i] = 'C'
            elif ilist[i] == 2: cl[i] = 'G'
            else:               cl[i] = 'U' # ilist[i] == 3

        cwrk = ''.join(ilist[:5])

        for i in xrange(100):
            if params.triloops[i] == '':
                break
            if params.triloops[i] == cwrk:
                el += params.dG_triloops[i]

    # Tetraloops
    #      1 2 3 4 5 6      #
    #  5'  A W X Y Z U   3' #
    #        L O O P        #

    elif nl == 4:

        for i in xrange(6):
            if   ilist[i] == 0: cl[i] = 'A'
            elif ilist[i] == 1: cl[i] = 'C'
            elif ilist[i] == 2: cl[i] = 'G'
            else:               cl[i] = 'U' # ilist[i] == 3

        cwrk = ''.join(ilist[:6])

        for i in xrange(100):
            if params.tetraloops[i] == '':
                break
            if params.tetraloops[i] == cwrk:
                el += params.dG_tetraloops[i]

    # Hexaloops
    #      0 1 2 3 4 5 6 7      #
    #  5'  A Q V W X Y Z U   3' #
    #        - L O O P -        #

    else: # nl == 6:

        for i in xrange(8):
            if   ilist[i] == 0: cl[i] = 'A'
            elif ilist[i] == 1: cl[i] = 'C'
            elif ilist[i] == 2: cl[i] = 'G'
            else:               cl[i] = 'U' # ilist[i] == 3

        cwrk = ''.join(ilist[:6])

        for i in xrange(100):
            if params.hexaloops[i] == '':
                break
            if params.hexaloops[i] == cwrk:
                el += params.dG_hexaloops[i]
    
    return el