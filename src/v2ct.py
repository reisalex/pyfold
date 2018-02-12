"""
Subroutine: V2CT (IBSP,FLD,JOB,N)

Description: Converts a secondary structure from either Vienna format
       to CT format, or CT to Vienna.

Arguments:

      IBSP - Base-pair information for the RNA.
       FLD - Vienna fold for the RNA.
       JOB - Conversion job to perform.
         JOB = 'V' Convert to Vienna Format
         JOB = 'C' Convert to CT Format.
       N - Number of nucleotides in the RNA.

History:
Version     Date            Comment
--------    -------         --------------------
      09/28/2017      Original Code

Dependencies:

Author(s): Alex Reis
       Copyright (c) 2017 (Please refer to LICENCE)
"""

def V2CT(fld):

  # INTEGERS
  # n,ibsp(n),i,j,count,istack(n)

  # STRING
  # fld(n)

  n = len(fld)

  s = 0

  ibsp   = [0]*n
  istack = [0]*n

  for i in range(0,n):

    if fld[i] == '.':
      ibsp[i] = 0

    elif fld[i] == '(':
      s += 1
      istack[s] = i

    elif fld[i] == ')':
      j = istack[s]
      ibsp[j] = i
      ibsp[i] = j
      s -= 1

  return ibsp

def CT2V(ibsp,n):

  # INTEGERS
  # n,ibsp(n),i,j,count,istack(n)

  # STRING
  # fld(n)

  for i in range(0,n):

    if ( ibsp[i] == 0 ):
      fld[i] = '.'
    else:
      j = ibsp[i]
      if ( j > i ):
        fld[i] = '('
        fld[j] = ')'
      else:
        fld[i] = ')'
        fld[j] = '('

  return fld