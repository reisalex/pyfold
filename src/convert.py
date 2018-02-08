"""
Subroutine: CONVERT (GSEQ,ISEQ,NN)

Description: Converts a character sequence of A C G U into numerical code.

Method: Converts a nucleic acid sequence to the following numerical code:

        A = 1
        C = 2
        G = 3
      T/U = 4

Arguments:

          GSEQ - (INPUT) Array of length NN containing the sequence
                 of characters to be converted into numbers.
          ISEQ - (OUTPUT) Array of length NN containing the sequence
                 in the numerical code.
            NN - Number of characters in the sequence.

History:
Version     Date            Comment
--------    -------         --------------------
            09/28/2017      Original Code

Dependencies:

Author(s): Alex Reis
           Copyright (c) 2017 (Please refer to LICENCE)
"""

def CONVERT(gseq,iseq,nn):

  # Arguments
  # nn, iseq(nn), gseq(nn)

  # Variables
  # i

  iseq = [0]*nn

  for i in range(0,nn):
    
    if   gseq[i] == 'A': iseq[i] = 1
    elif gseq[i] == 'C': iseq[i] = 2
    elif gseq[i] == 'G': iseq[i] = 3
    elif gseq[i] == 'U': iseq[i] = 4
    else: raise Exception('Bad character in RNA sequence.')

  return iseq