"""
Subroutine: RANDOM (ISEED)

Description: Generates a random number of the interval [0,1] given an
             integer seed.

Method: Efficient Fortran Programming, John Wiley and Sons, New York (1990)
        pp. 17-18 Library of Congress code QA76.73.F25 K78

Arguments:

    ISEED - Four byte integer seed

History:
Version     Date            Comment
--------    -------         --------------------
            09/28/2017      Original Code

Dependencies:

Author(s): Alex Reis
           Copyright (c) 2017 (Please refer to LICENCE)
"""

def RANDOM(iseed):

    # INTEGERS
    # iseed, hi,lo,test,a,m,q,r

    # FLOAT
    # random

    a = 16807
    m = 2147483647
    q = 127773
    r = 2836

    hi = int(iseed/q)
    lo = iseed % q

    test = a * lo - r * hi

    if ( test > 0 ):
        iseed = test
    else:
        iseed = test + m

    random = float(iseed) / float(m)

    return random,iseed