import os
import sys
import numpy as np

def num(s):
    try:
        return float(s)
    except ValueError:
        return 0.00

# nucleotide to index conversion dict
convert = {'A': 0,'C': 1, 'G': 2, 'U': 3}

# ordered list of possible base pairs used in the parameter files
# basepairs = ['CG','GC','GU','UG','AU','UA']
basepairs = [ (1,2), (2,1), (2,3), (3,2), (0,3), (3,0) ]

class FreeEnergyParameters:
    def __init__(self):
        self.dG_stack      = np.zeros(shape=tuple([4]*4), dtype=float)
        self.dG_stackh     = np.zeros(shape=tuple([4]*4), dtype=float)
        self.dG_stacki     = np.zeros(shape=tuple([4]*4), dtype=float)
        self.dG_dangle5    = np.zeros(shape=tuple([4]*3), dtype=float)
        self.dG_dangle3    = np.zeros(shape=tuple([4]*3), dtype=float)
        self.dG_int11      = np.zeros(shape=tuple([4]*6), dtype=float)
        self.dG_int21      = np.zeros(shape=tuple([4]*7), dtype=float)
        self.dG_int22      = np.zeros(shape=tuple([4]*8), dtype=float)
        self.dG_hloop      = np.zeros(30, dtype=float)
        self.dG_bulge      = np.zeros(30, dtype=float)
        self.dG_iloop      = np.zeros(30, dtype=float)
        self.dG_triloops   = np.zeros(100, dtype=float)
        self.dG_tetraloops = np.zeros(100, dtype=float)
        self.dG_hexaloops  = np.zeros(100, dtype=float)

        self.dH_stack      = np.zeros(shape=tuple([4]*4), dtype=float)
        self.dH_stackh     = np.zeros(shape=tuple([4]*4), dtype=float)
        self.dH_stacki     = np.zeros(shape=tuple([4]*4), dtype=float)
        self.dH_dangle5    = np.zeros(shape=tuple([4]*3), dtype=float)
        self.dH_dangle3    = np.zeros(shape=tuple([4]*3), dtype=float)
        self.dH_int11      = np.zeros(shape=tuple([4]*6), dtype=float)
        self.dH_int21      = np.zeros(shape=tuple([4]*7), dtype=float)
        self.dH_int22      = np.zeros(shape=tuple([4]*8), dtype=float)
        self.dH_hloop      = np.zeros(30, dtype=float)
        self.dH_bulge      = np.zeros(30, dtype=float)
        self.dH_iloop      = np.zeros(30, dtype=float)
        self.dH_triloops   = np.zeros(100, dtype=float)
        self.dH_tetraloops = np.zeros(100, dtype=float)
        self.dH_hexaloops  = np.zeros(100, dtype=float)

        self.triloops  = ['']*100
        self.tetraloops = ['']*100
        self.hexaloops = ['']*100

        # Other added parameters - see bottom of readpar()
        self.dG_bonuses    = np.zeros(6, dtype=float)
        self.dG_AU         = 0.0
        self.dG_asym       = np.zeros(6, dtype=float)

        self.dH_bonuses    = np.zeros(6, dtype=float)


def readpar(paramfile):

    params = FreeEnergyParameters()

    with open(paramfile,'r') as f:

        data = {}
        current = None
        for line in map(str.strip,f):
            if line.startswith('# '):
                data[line] = []
                current = line
            elif current:
                data[current].append(line)

    # ===============================================================
    # TSTACK
    # 5' (1) A X (3) 3'
    # 3' (2) U Y (4) 5'

    # *.par follows a 12(column):43(row) format instead of 12:34 used here
    
    readlines = iter(data['# stack'])
    next(readlines) # skip table column headers line

    # process stack free energies
    for i in xrange(6):
        line = next(readlines)
        values = map(num,line.split()[:6])
        for j in xrange(6):
            XY = basepairs[i]
            WZ = basepairs[j]
            params.dG_stack[ XY[0], XY[1], WZ[1], WZ[0] ] = values[j]

    params.dG_stack /= 100.0

    # skip lines before `# stack_enthalpies`
    readlines = iter(data['# stack_enthalpies'])
    next(readlines)

    for i in xrange(6):
        line = next(readlines)
        values = map(num,line.split()[:6])
        for j in xrange(6):
            XY = basepairs[i]
            WZ = basepairs[j]
            params.dH_stack[ XY[0], XY[1], WZ[1], WZ[0] ] = values[j]

    params.dH_stack /= 100.0

    # ===============================================================
    # TSTACKH
    # Stacking interaction of the 2 nt over a closing bp of a hairpin loop
    # 5' (1) A X (3) 3' LOOP!
    # 3' (2) U Y (4) 5' LOOP!

    readlines = iter(data['# mismatch_hairpin'])

    for i in xrange(6):
        for j in xrange(-1,4):
            line = next(readlines)
            if j == -1:
                continue
            values = map(num,line.split()[1:5])
            XY = basepairs[i]
            params.dG_stackh[ XY[0], XY[1],j, : ] = values

    params.dG_stackh /= 100.0

    readlines = iter(data['# mismatch_hairpin_enthalpies'])

    for i in xrange(6):
        for j in xrange(-1,4):
            line = next(readlines)
            if j == -1:
                continue
            values = map(num,line.split()[1:5])
            XY = basepairs[i]
            params.dH_stackh[ XY[0], XY[1],j, : ] = values

    params.dH_stackh /= 100.0

    # ===============================================================
    # TSTACKI
    # Stacking interaction of the 2 nt over a closing bp of a hairpin loop
    # 5' (1) A X (3) 3' INTERNAL LOOP!
    # 3' (2) U Y (4) 5' INTERNAL LOOP!

    readlines = iter(data['# mismatch_interior'])

    for i in xrange(6):
        for j in xrange(-1,4):
            line = next(readlines)
            if j == -1:
                continue
            values = map(num,line.split()[1:5])
            XY = basepairs[i]
            params.dG_stacki[ XY[0], XY[1],j, : ] = values

    params.dG_stacki /= 100.0

    readlines = iter(data['# mismatch_interior_enthalpies'])

    for i in xrange(6):
        for j in xrange(-1,4):
            line = next(readlines)
            if j == -1:
                continue
            values = map(num,line.split()[1:5])
            XY = basepairs[i]
            params.dH_stacki[ XY[0], XY[1],j, : ] = values

    params.dH_stacki /= 100.0

    # ===============================================================
    # TDANGLE5
    # 5' (1) A       3'
    # 3' (2) U X (3) 5'

    # *.par follows a 21(row):3(column) format compared to the 12:3 used here

    readlines = iter(data['# dangle5'])
    next(readlines)

    for i in xrange(6):
        line = next(readlines)
        values = map(num,line.split()[1:5])
        XY = basepairs[i]
        params.dG_dangle5[ XY[1], XY[0], : ] = values

    params.dG_dangle5 /= 100.0

    readlines = iter(data['# dangle5_enthalpies'])
    next(readlines)

    for i in xrange(6):
        line = next(readlines)
        values = map(num,line.split()[1:5])
        XY = basepairs[i]
        params.dH_dangle5[ XY[1], XY[0], : ] = values

    params.dH_dangle5 /= 100.0

    # ===============================================================
    # TDANGLE3
    # 5' (1) A X (3) 3'
    # 3' (2) U       5'

    # *.par follows a 21(row):3(column) format compared to the 12:3 used here

    readlines = iter(data['# dangle3'])
    next(readlines)

    for i in xrange(6):
        line = next(readlines)
        values = map(num,line.split()[1:5])
        XY = basepairs[i]
        params.dG_dangle3[ XY[1], XY[0], : ] = values

    params.dG_dangle3 /= 100.0

    readlines = iter(data['# dangle3_enthalpies'])
    next(readlines)

    for i in xrange(6):
        line = next(readlines)
        values = map(num,line.split()[1:5])
        XY = basepairs[i]
        params.dH_dangle3[ XY[1], XY[0], : ] = values

    params.dH_dangle3 /= 100.0

    # ========================================================================
    # TINT11
    #         (3)
    # 5' (1) A . X (5) 3'
    # 3' (2) U . Y (6) 5'
    #         (4)

    # *.par follows as 12:56:43 format instead of 12:34:56 as shown above
    # converting to 12:34:56 for now to conform to kfold format

    readlines = iter(data['# int11'])

    line = next(readlines)
    values = line.split()
    version = 1 if len(values) == 8 else 2

    if version == 1:
        while values:
            nts = values[-2].split(',')
            if nts[0] == 'NN' or nts[1] == 'NN' or nts[2] == 'N':
                pass
            else:
                values = values[1:5]
                X = [convert[n] for n in "".join(nts)]
                params.dG_int11[ X[0], X[1], X[4], :, X[3], X[2] ] = values
            line = next(readlines)
            values = line.split()
    else: # version == 2
        while values:
            nts = values[1].replace('..','')
            if len(nts) < 4:
                for _ in xrange(6):
                    line = next(readlines)
            else:
                X = [convert[n] for n in "".join(nts)]
                line = next(readlines) # skip N row
                for i in xrange(4):
                    line = next(readlines)
                    values = line.split()[1:]
                    params.dG_int11[ X[0], X[1], i, :, X[3], X[2] ] = values
                line = next(readlines)
            values = line.split()

    params.dG_int11 /= 100.0

    # ========================================================================
    # TINT12
    #         (3)
    # 5' (1) A .    X (6) 3'
    # 3' (2) U .  . Y (7) 5'
    #         (4)(5)

    # # int21
    # rna_turner2004.par table format:
    # ?     A     C     G     U
    # 230   230   230   110   230    /* CG,CG,A,C */
    # 5' (1) C A   G (7) 3'
    # 3' (2) G N C C (8) 5'    

    readlines = iter(data['# int21'])

    line = next(readlines)
    values = line.split()
    version = 1 if len(values) == 8 else 2

    if version == 1:
        while values:
            nts = values[-2].split(',')
            if nts[0] == 'NN' or nts[1] == 'NN' or nts[2] == 'N' or nts[3] == 'N':
                pass
            else:
                values = values[1:5]
                X = [convert[n] for n in "".join(nts)]
                params.dG_int21[ X[0], X[1], X[4], :, X[5], X[3], X[2] ] = values
            line = next(readlines)
            values = line.split()
    else: # version 2
        while values:
            nts = values[1].replace('.','')
            if len(nts) < 5 or nts.find('@') > -1:
                for _ in xrange(6):
                    line = next(readlines)
            else:
                X = [convert[n] for n in "".join(nts)]
                line = next(readlines) # skip N row
                for i in xrange(4):
                    line = next(readlines)
                    values = line.split()[1:]
                    params.dG_int21[ X[0], X[1], X[2], :, i, X[4], X[3] ] = values
                line = next(readlines)
            values = line.split()

    params.dG_int21 /= 100.0

    # ========================================================================
    # TINT22
    #         (3)(5)
    # 5' (1) A .  . X (7) 3'
    # 3' (2) U .  . Y (8) 5'
    #         (4)(6)

    # # int22
    # rna_turner2004.par table format:
    # A   C   G   U
    # 120 160 20  160 /* CG,CG,A,A,C */
    # 5' (1) C A A G (7) 3'
    # 3' (2) G N C C (8) 5'

    # turner1999 and andronescu2007 format:
    # /* CG.AA..CG */
    # 4x4 table (ACGU x ACGU)
    # 5' (1) C A A G (7) 3'
    # 3' (2) G N N C (8) 5'

    readlines = iter(data['# int22'])

    line = next(readlines)
    values = line.split()
    version = 1 if len(values) == 7 else 2

    if version == 1:
        while values:
            nts = values[-2].replace(',','')
            values = values[:4]
            X = [convert[n] for n in nts]
            params.dG_int22[ X[0], X[1], X[4], :, X[5], X[6], X[3], X[2] ] = values
            line = next(readlines)
            values = line.split()
    else: # version 2
        while values:
            nts = values[1].replace('.','')
            X = [convert[n] for n in "".join(nts)]
            for i in xrange(4):
                line = next(readlines)
                values = map(num,line.split())
                params.dG_int22[ X[0], X[1], X[2], :, X[3], i, X[5], X[4] ] = values
            line = next(readlines)
            values = line.split()

    params.dG_int22 /= 100.0

    # ========================================================================
    # Destabilizing energy of a HAIRPIN loop

    readlines = iter(data['# hairpin'])
    line = next(readlines)
    values = line.split()[1:]

    elh = []
    while values:
        elh.extend(map(num,values))
        line = next(readlines)
        values = line.split()
    params.dG_hloop[:] = elh
    params.dG_hloop /= 100.0


    readlines = iter(data['# hairpin_enthalpies'])
    line = next(readlines)
    values = line.split()[1:]

    hlh = []
    while values:
        hlh.extend(map(num,values))
        line = next(readlines)
        values = line.split()
    params.dH_hloop[:] = hlh
    params.dH_hloop /= 100.0

    # ========================================================================
    # Destabilizing energy of a BULGE loop

    readlines = iter(data['# bulge'])
    line = next(readlines)
    values = line.split()[1:]

    elb = []
    while values:
        elb.extend(map(num,values))
        line = next(readlines)
        values = line.split()
    params.dG_bulge[:] = elb
    params.dG_bulge /= 100.0


    readlines = iter(data['# bulge_enthalpies'])
    line = next(readlines)
    values = line.split()[1:]

    hlb = []
    while values:
        hlb.extend(map(num,values))
        line = next(readlines)
        values = line.split()
    params.dH_bulge[:] = hlb
    params.dH_bulge /= 100.0

    # ========================================================================
    # Destabilizing energy of a INTERIOR loop

    readlines = iter(data['# interior'])
    line = next(readlines)
    values = line.split()[1:]

    eli = []
    while values:
        eli.extend(map(num,values))
        line = next(readlines)
        values = line.split()
    params.dG_iloop[:] = eli
    params.dG_iloop /= 100.0


    readlines = iter(data['# interior_enthalpies'])
    line = next(readlines)
    values = line.split()[1:]

    hli = []
    while values:
        hli.extend(map(num,values))
        line = next(readlines)
        values = line.split()
    params.dH_iloop[:] = hli
    params.dH_iloop /= 100.0

    # ========================================================================
    # TLOOP
    # Tetra-loop bonus energies
    #    1 2 3 4 5 6
    # 5' A W X Y Z U 3'
    #      L O O P

    readlines = iter(data['# Triloops'])
    i = 0
    for line in readlines:
        if line == '':
            break
        values = line.split()
        params.triloops[i]   = values[0]
        params.dG_triloops[i] = num(values[1])
        params.dH_triloops[i] = num(values[2])
        i += 1

    readlines = iter(data['# Tetraloops'])
    i = 0
    for line in readlines:
        if line == '':
            break
        values = line.split()
        params.tetraloops[i]   = values[0]
        params.dG_tetraloops[i] = num(values[1])
        params.dH_tetraloops[i] = num(values[2])
        i += 1

    readlines = iter(data['# Hexaloops'])
    i = 0
    for line in readlines:
        if line == '':
            break
        values = line.split()
        params.hexaloops[i]   = values[0]
        params.dG_hexaloops[i] = num(values[1])
        params.dH_hexaloops[i] = num(values[2])
        i += 1

    # ========================================================================
    # Asymmetry of internal loop penalty (Ninio equation)

    raise Exception('Write import here')

    # ========================================================================
    # Hairpin loop bonuses and penalties

    if paramfile == 'rna_turner1999.par':

        params.dG_bonuses[0] = -0.8 # UU or GA first mismatch
        params.dG_bonuses[1] =  0.0 # GG first mismatch
        params.dG_bonuses[2] = -2.2 # special GU closure
        params.dG_bonuses[3] =  1.4 # C3 loop
        params.dG_bonuses[4] =  0.3 # All-C loop, A
        params.dG_bonuses[5] =  1.6 # All-C loop, B

        # values not provided on Turner website
        # params.dH_bonuses[0] = # UU or GA first mismatch
        # params.dH_bonuses[1] = # GG first mismatch
        # params.dH_bonuses[2] = # special GU closure
        # params.dH_bonuses[3] = # C3 loop
        # params.dH_bonuses[4] = # All-C loop, A
        # params.dH_bonuses[5] = # All-C loop, B

        params.dG_AU      = 0.5
        params.dG_asym[:] = 0.5
        params.dG_maxasym = 3.0

    elif paramfile == 'rna_turner20004.par':

        params.dG_bonuses[0] = -0.9 # UU or GA first mismatch
        params.dG_bonuses[1] = -0.8 # GG first mismatch
        params.dG_bonuses[2] = -2.2 # special GU closure
        params.dG_bonuses[3] =  1.5 # C3 loop
        params.dG_bonuses[4] =  0.3 # All-C loop, A
        params.dG_bonuses[5] =  1.6 # All-C loop, B

        params.dH_bonuses[0] = -5.8 # UU or GA first mismatch
        params.dH_bonuses[1] =  0.0 # GG first mismatch
        params.dH_bonuses[2] = -14.8 # special GU closure
        params.dH_bonuses[3] =  18.6 # C3 loop
        params.dH_bonuses[4] =  3.4 # All-C loop, A
        params.dH_bonuses[5] =  17.6 # All-C loop, B

        params.dG_AU      = 0.7
        params.dG_asym[:] = 0.6
        params.dG_maxasym = 3.0


    return params

if __name__ == "__main__":
    readpar(paramfile='rna_turner1999.par')
    readpar(paramfile='rna_turner2004.par')