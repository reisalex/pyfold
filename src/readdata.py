"""
Function: READDATA

Description: Reads in the RNA energy data files and sets up some tables.

History:
Version     Date            Comment
--------    -------         --------------------
            09/28/2017      Original Code

Dependencies:

Author(s): Alex Reis
           Copyright (c) 2017 (Please refer to LICENCE)
"""

from rnavar import eau,eaup,iwc

def READDATA():

    # VARIABLES

    # INTEGERS
    # i,j,k

    #=== A=0,C=1,G=2,U=3 ===#

    iwc = [ [0]*4 for _ in range(4) ]
    iwc[0][3] = 1
    iwc[3][0] = 1
    iwc[1][2] = 1
    iwc[2][1] = 1
    iwc[2][3] = 1
    iwc[3][2] = 1

    eaup = [ [0.0e0]*4 for _ in range(4) ]
    eaup[0][3] = eau
    eaup[3][0] = eau
    eaup[2][3] = eau
    eaup[3][2] = eau

    return
