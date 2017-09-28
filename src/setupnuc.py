"""
Subroutine: SETUPNUC

Description: Creates a table of nucleation probabilities between pairs of
             nucleotides based on a worm like chain model.

            (1) Toan et al. J. Phys. Chem. B 112, 6094-6106 (2008).
            (2) S. Kuznetsov and A. Ansari "A kinetic zipper model with
            interchain interactions applied to nucleic acid hairpin
            folding kinetics", Biophys J 102, 1001-111 (2012).

History:
Version     Date            Comment
--------    -------         --------------------
            09/28/2017      Original Code

Dependencies:

Author(s): Alex Reis
           Copyright (c) 2017 (Please refer to LICENCE)
"""
import math
from rnavar import pnuc,beta

def SETUPNUC(n):

    # Variables
    # INTEGER
    # n,i
    
    # FLOAT
    # x,xi,e,c,c2,xp

    c = 0.1785714290e0
    c2 = 3.9274668195e2 # rate in (1/uS)
    xp = 4.0000000000e0

    pnuc[:] = 0.0e0

    for i in range(4,n):

        x = c * float(i-1)
        xi = 1.0e0 / x

        if ( x <= xp ):

            e = -7.0270e0 * xi + 0.4920e0 * x
            x = 84.90e0 * ( xi ** 5.50e0 )
            x = c2 * c2 / beta

            pnuc[i] = x * math.exp(e)

        else:

            x = xi ** 2.0

            e = 1.0e0 - 0.6250e0 * xi - 0.12343750e0 * x
            x = c2 * x / beta

            pnuc[i] = x * e


    return pnuc
    

