"""
Module: DELTAG_HD (RNA,II,JJ,KK)

Description: Computes the difference in free energy of an RNA loop due
             to a nucleotide diffusion along the helix i.e. the ii-jj
             base pair will shift to either ii-kk or kk-jj. The energy
             change is calculated using the empirical INN model.

Arguments:

          RNA - Class containing information about the RNA fold and
                the RNA sequence.
           II - Nucleotide position of the 5' most nucleotide.
           JJ - Nucleotide position of the 3' most nucleotide.
           KK - Nucleotide position of the single-stranded nucleotiode
                that either ii/jj in the base-pair ii-jj will swap with.

History:
Version     Date            Comment
--------    -------         --------------------
            09/28/2017      Original Code

Dependencies:

Author(s): Alex Reis
           Copyright (c) 2017 (Please refer to LICENCE)
"""
