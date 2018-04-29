"""
Module: RNAVAR

Description: Contains the variables, work arrays, and parameters needed for computing
             the kinetics of folding for an RNA molecule. 

History:
Version     Date            Comment
--------    -------         --------------------
            09/28/2017      Original Code

Dependencies:

Author(s): Alex Reis
           Copyright (c) 2017 (Please refer to LICENCE)
"""

#Global
# global pnuc,iwc,eaup
# global beta
# global gcons,temp
# global rateh,ratem,rated
# global mxnt
# global em,eh,es,eau
# global parameters

# Allocatable arrays
pnuc = []

# Variables
beta = 0.16225023135094183147e1

iwc = [] # 4x4 matrix integer
eaup = [] # 4x4 matrix real

# Parameters
gcons = 1.987206e-3
temp = 310.15e0

rateh = 1.0e2 # helix extension (1/uS)
ratem = 5.0e0 # helix morphing  (1/uS)
rated = 5.0e0 # helix diffusion (1/uS)

mxnt = 10000

em = 10.10e0
eh = -0.30e0
es = -0.30e0
eau = 0.50e0


class EnergyParameters(object):

    # available_files = [   "dna_mathews1999.par",
    #                        "rna_turner1999.par",
    #                       "dna_mathews2004.par",
    #                        "rna_turner2004.par",
    #                    "rna_andronescu2007.par"]

    temp = 310.15 # default temperature (37.0 C)

    #=== ViennaRNA .par files ===#

    '''
    dG_stack_37C
    dH_stack

    dG_mismatch_hairpin_37C
    dH_mismatch_hairpin

    dG_mismatch_interior_37C
    dH_mismatch_interior

    dG_mismatch_interior_1n_37C
    dH_mismatch_interior_1n

    dG_mismatch_interior_23_37C
    dH_mismatch_interior_23

    dG_mismatch_multi_37C
    dH_mismatch_multi

    dG_mismatch_exterior_37C
    dH_mismatch_exterior

    dG_dangle5_37C
    dH_dangle5

    dG_dangle3_37C
    dH_dangle3

    dG_int11_37C
    dH_int11

    dG_int21_37C
    dH_int21

    dG_int22_37C
    dH_int22

    hairpin
    hairpin_enthalpies

    bulge
    bulge_enthalpies

    interior
    interior_enthalpies

    ML_params
    NINIO

    Misc
    Hexaloops

    Tetraloops
    Triloops
    '''


    #=== KFOLD REQUIRED ENERGIES ===#

    '''
    ========================================================================
    EBULGE
    Compute the energy of an RNA bulge with 2 helices
    EB = E_entropic + E_stack + E_stack + E_asymmetry
    5' (I) X ... W (IP) 3'
    3' (J) Y ... Z (JP) 5'
    elb,eli (entropic term)

    EHAIR
    Compute the energy of an RNA hairpin turn
    EH = E_entropic + E_stack + E_bonus + E_penalty
    5' (I) X ... loop 3'
    3' (J) Y ... loop 5'
    elh (entropic term)

    ========================================================================
    TSTACK
    5' (1) A X (3) 3'
    3' (2) U Y (4) 5'

    ========================================================================
    TDANGLE3
    5' (1) A X (3) 3'
    3' (2) U       5'

    ========================================================================
    TDANGLE5
    5' (1) A       3'
    3' (2) U X (3) 5'

    ========================================================================
    TINT11
            (3)
    5' (1) A . X (5) 3'
    3' (2) U . Y (6) 5'
            (4)

    ========================================================================
    TINT12
            (3)
    5' (1) A .    X (6) 3'
    3' (2) U .  . Y (7) 5'
            (4)(5)

    # int21
    rna_turner2004.par table format:
    ?     A     C     G     U
    230   230   230   110   230    /* CG,CG,A,C */
    5' (1) C A   G (7) 3'
    3' (2) G N C C (8) 5'    

    ========================================================================
    TINT22
            (3)(5)
    5' (1) A .  . X (7) 3'
    3' (2) U .  . Y (8) 5'
            (4)(6)

    # int22
    rna_turner2004.par table format:
    A   C   G   U
    120 160 20  160 /* CG,CG,A,A,C */
    5' (1) C A A G (7) 3'
    3' (2) G N C C (8) 5'

    turner1999 and andronescu2007 format:
    /* CG.AA..CG */
    4x4 table (ACGU x ACGU)
    5' (1) C A A G (7) 3'
    3' (2) G N N C (8) 5'

    ========================================================================
    TLOOP
    Tetra-loop bonus energies
       1 2 3 4 5 6
    5' A W X Y Z U 3'
         L O O P

    # Tetraloops
    sequence   dG      dH
    CAACGG     550     690

    # Add hexa- and tri- loops?

    ========================================================================
    TSTACK
    Stacking interaction of the 2 nt over a closing bp of a hairpin loop
    5' (1) A X (3) 3' LOOP! <--- Double check loop position
    3' (2) U Y (4) 5' LOOP!
    
    # stack

    '''


    def __init__(self):
        pass

    def set_temp(self,temperature):
        self.temp = temperature + 273.15 # convert to Kelvin


