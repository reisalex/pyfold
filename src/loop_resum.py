"""
Subroutine: LOOP_RESUM (RNA,INDX)

Description: Resums the partial sum table of transition rates and total flux
       when the transition rate for a single loop element (#indx)
       has changed.

Method: Recomputes the total flux in LOG_2(N) time using a partial sum table.

Arguments:
    
       R - Class structure containing information on the
         RNA secondary structure and possible reactions.
      INDX - The indx number of the loop element

History:
Version     Date            Comment
--------    -------         --------------------
      09/28/2017      Original Code

Dependencies:

Author(s): Alex Reis
       Copyright (c) 2017 (Please refer to LICENCE)
"""

def LOOP_RESUM(rna,indx):

  # INTEGER
  # i,j,k,nsum
  # n,n1,n2

  nsum = rna.nsum

  #=== Resum Partial Sum Table ===#

  n = 1
  n1 = 2
  n2 = 4

  if ( indx % 2 == 1 ): i = indx
  if ( indx % 2 == 0 ): i = indx - 1

  rna.psum[i] = rna.ptot[i] + rna.ptot[i+1]

  while ( n1 < nsum ):

    i = int(i/n2) * n2 + n1

    j = i - n
    k = i + n

    rna.psum[i] = rna.psum[j] + rna.psum[k]

    n = n1
    n1 = n2
    n2 = 2 * n2

  #endwhile

  return rna
