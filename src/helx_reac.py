"""
Subroutine: LOOP_FIRE (RNA,INDX,AMAX)

Description: Finds the reaction J in the loop such that S >= AMAX where
             S is the partial sum of the reaction rates for reactions
             {1,J} in the loop. AMAX is determined from the SSA protocol.
             Once reaction J is found, this reaction is "fired" and the
             reactions for neighboring loop elements updated.

Arguments:
        
             R - Class structure containing information on the
                 RNA secondary structure and possible reactions.
          INDX - The indx number of the loop element that a reaction
                 will be choosen.
          AMAX - The reaction "threshold" determined from the SSA protocol.

History:
Version     Date            Comment
--------    -------         --------------------
            09/28/2017      Original Code

Dependencies:

Author(s): Alex Reis
           Copyright (c) 2017 (Please refer to LICENCE)
"""