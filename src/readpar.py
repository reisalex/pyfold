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
        self.dG_stack   = np.zeros(shape=(4,4,4,4),         dtype=float)
        self.dH_stack   = np.zeros(shape=(4,4,4,4),         dtype=float)
        self.dG_dangle5 = np.zeros(shape=(4,4,4),           dtype=float)
        self.dH_dangle5 = np.zeros(shape=(4,4,4),           dtype=float)
        self.dG_dangle3 = np.zeros(shape=(4,4,4),           dtype=float)
        self.dH_dangle3 = np.zeros(shape=(4,4,4),           dtype=float)
        self.dG_int11   = np.zeros(shape=(4,4,4,4,4,4),     dtype=float)
        self.dH_int11   = np.zeros(shape=(4,4,4,4,4,4),     dtype=float)
        self.dG_int21   = np.zeros(shape=(4,4,4,4,4,4,4),   dtype=float)
        self.dH_int21   = np.zeros(shape=(4,4,4,4,4,4,4),   dtype=float)
        self.dG_int22   = np.zeros(shape=(4,4,4,4,4,4,4,4), dtype=float)
        self.dH_int22   = np.zeros(shape=(4,4,4,4,4,4,4,4), dtype=float)


def readpar(paramfile):

    params = FreeEnergyParameters()

    with open(paramfile,'r') as f:
        readlines = iter(map(str.strip,f))
        current = None

        # ===============================================================
        # TSTACK
        # 5' (1) A X (3) 3'
        # 3' (2) U Y (4) 5'

        # *.par follows a 12(column):43(row) format instead of 12:34 used here

        # skip lines before `# stack`
        while current != '# stack':
            current = next(readlines)
            continue
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
        while current != '# stack_enthalpies':
            current = next(readlines)
            continue
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
        # TDANGLE5
        # 5' (1) A       3'
        # 3' (2) U X (3) 5'

        # *.par follows a 21(row):3(column) format compared to the 12:3 used here

        # skip lines before `# dangle5`
        while current != '# dangle5':
            current = next(readlines)
            continue
        next(readlines)

        for i in xrange(6):
            line = next(readlines)
            values = map(num,line.split()[1:5])
            XY = basepairs[i]
            params.dG_dangle5[ XY[1], XY[0], : ] = values

        params.dG_dangle5 /= 100.0

        # skip lines before `dangle5_enthalpies`
        while current != '# dangle5_enthalpies':
            current = next(readlines)
            continue
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

        # skip lines before `# dangle5`
        while current != '# dangle3':
            current = next(readlines)
            continue
        next(readlines)

        for i in xrange(6):
            line = next(readlines)
            values = map(num,line.split()[1:5])
            XY = basepairs[i]
            params.dG_dangle3[ XY[1], XY[0], : ] = values

        params.dG_dangle3 /= 100.0

        # skip lines before `dangle3_enthalpies`
        while current != '# dangle3_enthalpies':
            current = next(readlines)
            continue
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

        # skip lines before `# int11`
        while current != '# int11':
            current = next(readlines)
            continue

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

        # skip lines before `# int11`
        while current != '# int21':
            current = next(readlines)
            continue

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

        # skip lines before `# int11`
        while current != '# int22':
            current = next(readlines)
            continue

        line = next(readlines)
        values = line.split()
        version = 1 if len(values) == 7 else 2

        if version == 1:
            while values:
                line = next(readlines)
                values = line.split()
                # print values

            # while 1:
            #     nts = values[-2].split(',')
            #     if nts[0] == 'NN' or nts[1] == 'NN' or nts[2] == 'N' or nts[3] == 'N':
            #         pass
            #     else:
            #         values = values[1:5]
            #         X = [convert[n] for n in "".join(nts)]
            #         params.dG_int21[ X[0], X[1], X[4], :, X[5], X[3], X[2] ] = values
            #     line = next(readlines)
            #     values = line.split()
            #     if not values:
            #         break

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

        print params.dG_int22[0,3,0,0,0,1,1,2]
        # AAAC
        # UACG

        print params.dG_int22[0,3,0,0,1,1,1,2]
        # AACC
        # UACG

        print params.dG_int22[0,3,0,2,2,2,1,2]
        # AAGC
        # UGGG

        print params.dG_int22[0,3,3,3,3,2,1,2]
        # AUUC
        # UUGG

        memory = 0.0
        for attr in dir(params):
            if not attr.startswith('__'):
                memory += getattr(params,attr).nbytes
        print memory


        # ========================================================================
        # TLOOP
        # Tetra-loop bonus energies
        #    1 2 3 4 5 6
        # 5' A W X Y Z U 3'
        #      L O O P

        # # Tetraloops
        # sequence   dG      dH
        # CAACGG     550     690

        # # Add hexa- and tri- loops?

        # ========================================================================
        # TSTACK
        # Stacking interaction of the 2 nt over a closing bp of a hairpin loop
        # 5' (1) A X (3) 3' LOOP! <--- Double check loop position
        # 3' (2) U Y (4) 5' LOOP!
        
        # # stack

        # WRITE CODE TO ADJUST rna_turner1999.par BC `DEF` VALUES PRESENT

    # return params

if __name__ == "__main__":
    # readpar(paramfile='rna_turner1999.par')
    readpar(paramfile='rna_turner2004.par')