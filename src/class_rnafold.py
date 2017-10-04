"""
============================================================================
Module: CLASS_RNAFOLD

Description: A class structure containing subroutines and data elements
             required for computing transition probabilities between
             different RNA secondary structures.

History:
Version     Date            Comment
--------    -------         --------------------
            09/28/2017      Original Code

Dependencies:

Author(s): Alex Reis
           Copyright (c) 2017 (Please refer to LICENCE)
============================================================================           
"""


from rnavar import *

import loop_init
import loop_resum
import helx_reac
import loop_reac
import loop_fire

class RNA_STRUC(object):

    seq = "N"*mxnt

    iseq = [0]*mxnt
    ibsp = [0]*mxnt
    link = [0]*mxnt

    loop = [0]*mxnt
    nhlx = [0]*mxnt
    nsgl = [0]*mxnt

    n = 0
    nl = 0
    nsum = 0

    psum = [0.0]*mxnt
    ptot = [0.0]*mxnt

    wrk1 = [0.0]*mxnt
    wrk2 = [0.0]*mxnt

    def clear_loops(self):

        self.link = [0]*mxnt
        self.loop = [0]*mxnt
        self.nhlx = [0]*mxnt
        self.nsgl = [0]*mxnt

        self.wrk1 = [0.0]*mxnt
        self.wrk2 = [0.0]*mxnt

        self.psum = [0.0]*mxnt
        self.ptot = [0.0]*mxnt


    def LOOP_INIT(self):
        loop_init.LOOP_INIT(self)

    def LOOP_RESUM(self):
        loop_resum.LOOP_RESUM(self)

    def HELX_REAC(self):
        helx_reac.HELX_REAC(self)

    def LOOP_REAC(self,indx):
        loop_reac.LOOP_REAC(self,indx)

    def LOOP_FIRE(self,indx,amax):
        loop_fire.LOOP_FIRE(self,indx,amax)