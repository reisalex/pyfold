"""
Module: RNAVAR

Description: Contains the variables, work arrays, and parameters needed for computing
			 the kinetics of folding for an RNA molecule. 

History:
Version		Date			Comment
--------    -------     	--------------------
			09/28/2017		Original Code

Dependencies:

Author(s): Alex Reis
		   Copyright (c) 2017 (Please refer to LICENCE)
"""

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