"""
Program: PYFOLD

Description: A program for computing the folding kinetics of an RNA
             sequence using Turner free energies.

History:
Version     Date            Comment
--------    -------         --------------------
            09/28/2017      Original Code

Dependencies:

Author(s): Alex Reis
           Copyright (c) 2017 (Please refer to LICENCE)
"""

import re
from rnavar import mxnt
from ssareaction import SSAREACTION
from readdata import READDATA

def PYFOLD(seq,fld_start=None,fld_stop=None,nsim=1,tmax=1.0):

    # VARIABLES
    
    # STRINGS
    # seq,fld

    # ARRAYS
    # iseq,ibpi,ibpf

    # INTEGERS
    # i,j,k,n,nn,is,io,isim,nsim,narg,iseed

    # FLOAT
    # random,tstart,time,tout,tmax,dt

    # LOGICAL
    # istart,istop

    # DEFAULT SETTINGS

    READDATA()

    rna = RNA_STRUC()

    tstart = 0.0
    iseed = 61928712

    # Initial and final structure
    iseq = [0]*mxnt
    ibpi = [0]*mxnt
    ibpf = [0]*mxnt

    istart = False
    istop  = False

    nn = len(seq)

    # Check seq
    assert nn <= mxnt, "Error: Maximum number of nt = {}".format(mxnt)
    assert isinstance(seq,str), "Error: seq must be a string."
    seq = seq.upper().replace("T","U")
    exp = re.compile('[AGCU]',re.IGNORECASE)
    if exp.match(seq) == None:
        raise ValueError("Invalid letters found in sequence: {}. Only AGCU accepted.".format(seq))

    # Check and process start and stop structures
    exp = re.compile('[\.\(\)]')
    if not fld_start is None:
        assert isinstance(fld_start,str)
        if exp.match(fld_start) == None:
            raise ValueError("Invalid letters found in start structure: {}. Only .() accepted.".format(fld_start))
        istart = True
        ibpi = V2CT(fld_start,'C',nn)

    if not fld_stop is None:
        assert isinstance(fld_stop,str)
        if exp.match(fld_stop) == None:
            raise ValueError("Invalid letters found in stop structure: {}. Only .() accepted.".format(fld_stop))
        istop = True
        ibpf = V2CT(fld_stop,'C',nn)

    # Set up RNA
    SETUPNUC(nn) #!

    rna.seq = seq
    rna.iseq = CONVERT(seq,nn)
    rna.n = nn

    # Simulate RNA kinetics

    for isim in range(1,nsim+1):

        io = 1
        dt = 1.0e-2

        tout = dt
        time = tstart

        rna.ibsp = ibpi

        rna.LOOP_INIT()

        # Stochastic simulation
        while time < tmax:

            rna,iseed = SSAREACTION(rna,iseed,time,tout)

            # Increment tout
            if time > tout:
                tout = tout + dt
                io = io + 1
                if io > 9:
                    io = 1
                    dt = dt*10.0

            # Check for stop structure
            if istop:
                if ibpf == rna.ibsp:

                    # Write or whatever
                    return 

    return


