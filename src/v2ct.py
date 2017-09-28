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

def V2CT(ibsp=[],fld=[],job,n):

    # INTEGERS
    # n,ibsp(n),i,j,count,istack(n)

    # STRING
    # job,fld(n)

    #=== Convert to Vienna ===#
    if ( job == "V" ):

        assert (not ibsp is None), "Make sure to provide ibsp for job 'V'."

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



    #=== Convert to CT ===#
    if ( job == "C" ):

        assert (not fld is None), "Make sure to provide fld for job 'C'."

        count = 0
        ibsp[:] = 0
        istack[:] = 0

        for i in range(0,n):

            if fld[i] == '.':
                ibsp[i] = 0

            elif fld[i] == '(':
                count += 1
                istack[count] = i

            elif fld[i] == ')':
                j = istack[count]
                ibsp[j] = i
                ibsp[i] = j
                count -= 1

    return ibsp,fld