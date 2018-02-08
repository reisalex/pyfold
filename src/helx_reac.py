"""
Subroutine: HELX_REAC (RNA,INDX)

Description: Recomputes the reaction rates for opening one of the internal
             base-pairs within a continuous section of hliex.

Method:

              A-U           A-U
              x-x   --->   x   x
              U-A           U-A
Arguments:

            RNA  - Class structure containing information on the
                   RNA secondary structure and possible reactions.
            INDX - The indx number of one of the nucelotides which
                   belongs to the helix.

History:
Version     Date            Comment
--------    -------         --------------------
            09/28/2017      Original Code

Dependencies:

Author(s): Alex Reis
           Copyright (c) 2017 (Please refer to LICENCE)
"""

import math

def HELX_REAC(rna,indx):

  # Variables

  # INTEGER
  # ke,ip,jp,indx

  # FLOAT
  # x,dg,atot

  jp = indx

  ip = rna.ibsp[jp]
  jp = min(ip,jp)

  while ( rna.link[jp] == 0 ):
      jp += 1
  #endwhile

  indx = rna.link[jp]

  #=== Open BP Inside Helix ===#

  ke = rna.ibsp[jp]
  ip = rna.ibsp[jp]

  atot = rna.wrk1[jp]

  if ( rna.link[ke] == 0 ):

    ip += 1
    jp -= 1

    while ( rna.link[ip] == 0 ):

      dg = DELTAG_HI(rna,ip,jp)

      dg = dg / 2.0e0

      x = beta * dg
      x = math.exp(-x) * rateh

      rna.wrk1[ip] = x
      rna.wrk1[jp] = 0.0e0

      rna.wrk2[ip] = 0.0e0
      rna.wrk2[jp] = 0.0e0

      atot += x
      
      ip += 1
      ip -= 1

    #endwhile
  #endif

  rna.ptot[indx] = atot

  rna.LOOP_RESUM(indx)

  return

