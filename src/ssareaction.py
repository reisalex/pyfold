"""
Subroutine: SSAREACTION (RNA,ISEED,TIME,TOUT)

Description: Calculates an RNA folding reaction to fire based on the
             (S)tochastic (S)imulation (A)lgorithm of Gillespie.

See Gillespie, Daniel T. (1977). "Exact Stochastic Simulation
of Coupled Chemical Reactions". J. Phys. Chem. 81, 2340.

Arguments:

    RNA     - Class structure containing information on the
                RNA secondary structure and possible reactions.
    ISEED   - Integer seed for the random number generator.
    TIME    - Current Time
    TOUT    - Time to write next trajectory output.

History:
Version     Date            Comment
--------    -------         --------------------
            09/28/2017      Original Code

Dependencies:

Author(s): Alex Reis
           Copyright (c) 2017 (Please refer to LICENCE)
"""

import math

from rnavar import mxnt
from random import RANDOM


def SSAREACTION(rna,iseed,time,tout):

    # VARIABLES

    # REAL
    # e

    # INTEGER
    # i,j,k,n
    # n1,n2,indx

    # STRING
    # fld(mxnt)

    # FLOAT
    # r,tau,random
    # atot,amax


    #=== Total transition rate ===#
    n = rna.nsum / 2
    atot = rna.psum(n)

    #=== Compute time increment ===#
    r,iseed = RANDOM(iseed)

    tau = math.log(1.0/r)
    tau = tau / atot

    #=== Output current structure? ===#
    if time > tout:
        # do write things

    #=== Fire reaction ===#
    r,iseed = RANDOM(iseed)
    amax = r * atot

    #=== Find reaction to fire ===#
    i = n

    while (n % 2 == 0):
        n = n / 2
        j = i - n
        r = rna.psum(j)
        if (r >= amax):
            i = j
        else:
            i = i + n
            amax = amax - r

    #=== Choose between i and i+1 ===#
    r = rna.ptot(i)
    if ( r >= amax ):
        indx = i
    else:
        indx = i + 1
        amax = amax - r

    #=== Fire reaction ===#
    rna.LOOP_FIRE(indx,amax)

    return rna,iseed

