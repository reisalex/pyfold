"""
Subroutine: LOOP_REAC (RNA,INDX)

Description: Recomputes the possible reactions within a loop element and
       their corresponding rates.

! Method:  There are 6 possible reactions around the loop:
!
!          (1) Nucleation                   A
!                                          A A
!              A A A A A A U A A  -->  A A A-U A A
!
!          (2) Helix Extension
!                                        x-x
!              A A A x-x U A A  -->  A A A-U A A
!
!          (3) Helix Retraction
!                    x-x
!                A A A-U A A   -->  A A A x-x U A A
!
!          (4) Helix Morphing               x-x
!                  x-x  x-x                 U-A
!                A A-U  U-A A  -->  A A x-x U-A
!
!          (5) Defect Diffusion         x-x
!                    x-x                  U
!                A A A-U U A A -->  A A A-U A A
!
!          (6) Helix Open
!                         U-A       U-A
!                         x-x      x   x
!                         A-U  -->  A-U

Arguments:
    
       R - Class structure containing information on the
         RNA secondary structure and possible reactions.
      INDX - The indx number of the loop element that reactions
         will be calculated for.

History:
Version     Date            Comment
--------    -------         --------------------
      09/28/2017      Original Code

Dependencies:

Author(s): Alex Reis
       Copyright (c) 2017 (Please refer to LICENCE)
"""

import math

def LOOP_REAC(rna,indx):

  # VARIABLES

  # INTEGERS
  # indx,i,j,k,n,ip,jp,kp,hs,js
  # ks,ke,l,lmx,icnt,iloop
  # nt,nh,ns,mt,mh,ms,icase

  # FLOAT
  # x,dg,atot,rate

  i = rna.loop[indx]
  j = rna.ibsp[i]
  n = rna.n

  if ( i == n ): j = 1

  nh = rna.nhlx[indx]
  ns = rna.nsgl[indx]

  nt = ns + 2 * nh

  #--- Internal Loop iloop = 1 ---#
  #--- External Loop iloop = 0 ---#

  if ( i < j ): iloop = 1
  if ( i > j ): iloop = 0

  if ( iloop == 1 ):
    ks = i + 1
    ke = j
  else:
    ks = j
    ke = i
  #endif

  #=== COMPUTE REACTIONS ===#

  atot = 0.0e0

  k = ks
  icnt = 0

  while ( k <= ke ):

    #=== Nucleation Events ===#

    if ( rna.ibsp[k] == 0 ):

      rna.wrk1[k] = 0.0e0
      rna.wrk2[k] = 0.0e0

      x = 0.0e0

      l = 2
      lmx = nt / 2 + 1

      if ( nt % 2 == 0 ) and ( icnt+1 > lmx-1 ):
        lmx -= 1
      #endif

      if ( iloop == 0 ): lmx = nt - icnt

      kp = k + 1
      hs = rna.iseq[k]

      while ( l <= lmx ):

        if ( rna.ibsp[kp] == 0 ):

          js = rna.iseq[kp]

          if ( l > 4 and iwc[hs][js] == 1 ):
            x += pnuc[l]
          #endif

        else:

          l += 1
          kp = rna.ibsp[kp]
        #endif
        
        l += 1
        kp += 1

      #endwhile

      rna.wrk1[k] = x
      atot += x

    #endif

    #=== Helix Events ===#

    if ( rna.ibsp[k] > 0 ):

      ip = k
      jp = rna.ibsp[k]

      rna.wrk1[jp] = 0.0e0
      rna.wrk2[ip] = 0.0e0

      if ( rna.link[ip] == 0 ):
        rna.wrk1[ip] = 0.0e0
        rna.wrk2[jp] = 0.0e0
      #endif

      icase = 0

      #=== Helix Extension ===#

      if  ( ip > 1 and jp < n ) and \
        ( nh > 1 or ns > 4 ) and \
        ( rna.ibsp[ip-1] == 0 and rna.ibsp[jp+1] == 0 ):

        hs = rna.iseq[ip-1]
        js = rna.iseq[jp+1]

        if ( iwc[hs][js] == 1 ):

          icase = 1

          if ( iloop == 1 ):
            if ( nh == 2 and ns == 2 ):
              icase = 2
              if ( k == ke ): icase = 0
            elif ( k == ke ):
              icase = 3
            #endif
          #endif
        #endif
      #endif

      if ( icase > 0 ):

        dg = DELTAG_HE(rna,ip,jp,dg)

        dg = dg / 2.0e0

        x = beta * dg
        x = math.exp(-x) * rateh

        rna.wrk2[ip] = x
        atot += x

      #endif

      icase = 0

      #=== Helix Retraction ===#

      if ( ip != n and jp != 1 ):

        if ( rna.ibsp[ip+1] == jp-1 ):

          icase = 1

          if ( iloop == 1 ) and ( k == ke ):
            icase = 2
          #endif

          rate = rateh

        elif ( iloop == 0 or k != ke ):

          icase = 3

          l = rna.link[ip]
          mh = rna.nhlx[l]
          ms = rna.nsgl[l]

          mt = ms + 2 * mh

          if ( iloop == 1 ):
            l = min(nt,mt)
          else:
            l = mt
          #endif

          rate = pnuc[l]

        #endif
      #endif

      if ( icase > 0 ):

        dg = DELTAG_HR(rna,ip,jp,dg)

        if ( icase != 3 ):
          dg = dg / 2.0e0
        #endif

        x = beta * dg
        x = mat.exp(-x) * rate

        if ( icase == 2 ):
          rna.wrk1[ip] = x
        else:
          rna.wrk1[jp] = x
        #endif

        atot += x

      #endif

      icase = 0

      #=== Helix Morphing ===#

      if  ( iloop == 0 or nh > 2 ) and \
        ( ip > 1 and jp < n ):

        hs = rna.iseq[ip-1]
        js = rna.iseq[jp+1]

        if ( iwc[hs][js] == 1 ): icase = 1

        hs = rna.ibsp[ip-1]
        js = rna.ibsp[jp+1]

        if ( hs != 0 ) and ( rna.link[hs] != 0 ):
          icase = 0
        #endif

        if ( js != 0 ) and (rna.link[jp+1] != 0 ):
          icase = 0
        #endif

        if ( hs == 0 and js == 0 ):
          icase = 0
        #endif

      #endif

      if ( icase > 0 ):

        dg = DELTAG_HM(rna,ip,jp,dg)

        dg = dg / 2.0e0

        x = beta * dg
        x = math.exp(-x) * ratem

        atot += x

      #endif

      #=== Defect Diffusion ===#

      #--- PUSH ---#

      icase = 0

      if ( rna.link[ip] == 0 ):

        icase = 2

        if ( iloop == 1 ):
          if ( nh == 2 and ns == 1 ): icase = 3
          if ( nh == 1 and ns == 3 ): icase = 0
        #endif

      elif ( iloop == 0 or k != ke ):

        icase = 1

        if ( iloop == 1 ) and ( nh == 2 and ns == 1 ):
          icase = 4
        #endif

      #endif

      if ( icase > 0 ):

        #--- Push 5' end ---#

        kp = ip - 1

        if ( kp >= 1 ) and ( rna.ibsp[kp] == 0 ):

          hs = rna.iseq[kp]
          js = rna.iseq[jp]

          if ( iwc[hs][js] == 1 ):

            dg = DELTAG_HD(rna,ip,jp,kp)

            dg = dg / 2.0e0

            x = beta * dg
            x = math.exp(-x) * rated

            atot += x

          #endif
        #endif

        #--- Push 3' end ---#

        kp = jp + 1

        if ( kp <= n ) and ( rna.ibsp[kp] == 0 ):

          hs = rna.iseq[ip]
          js = rna.iseq[kp]

          if ( iwc[hs][js] == 1 ):

            dg = DELTAG_HD(rna,ip,jp,kp)

            dg = dg / 2.0e0

            x = beta * dg
            x = math.exp(-x) * rated

            atot += x

          #endif
        #endif
      #endif

      #--- PULL ---#

      icase = 0

      if  ( rna.link[ip] != 0 ) and \
        ( iloop == 0 or k != 0 ):

        l  = rna.link[ip]
        mh = rna.nhlx[l]
        ms = rna.nsgl[l]

        icase = 1

        if ( mh == 1 and ms == 3 ): icase = 0
        if ( mh == 2 and ms == 1 ): icase = 4

      #endif

      if ( icase > 0):

        #--- Pull 5' end ---#

        kp = ip + 1

        if ( rna.ibsp[kp] == 0 ):

          hs = rna.iseq[kp]
          js = rna.iseq[jp]

          if ( iwc[hs][js] == 1 ):

            dg = DELTAG_HD(rna,ip,jp,kp)

            dg = dg / 2.0e0

            x = beta * dg
            x = math.exp(-x) * rated

            atot += x

          #endif
        #endif

        #--- Pull 3' end ---#

        kp = jp - 1

        if ( rna.ibsp[kp] == 0 ):

          hs = rna.iseq[ip]
          js = rna.iseq[kp]

          if ( iwc[hs][js] == 1 ):

            dg = DELTAG_HD(rna,ip,jp,kp)

            dg = dg / 2.0e0

            x = beta * dg
            x = math.exp(-x) * rated

            atot += x

          #endif
        #endif
      #endif

      #=== Save Loop Reaction Rate ===#

      if ( iloop == 1 ): rna.wrk1[i] = atot

      #=== Open Internal Helix BP ===#

      if  ( rna.link[ip] == 0 ) and \
        ( iloop == 1 and k == ke ):

        hs = ip + 1
        js = jp - 1

        while ( rna.link[hs] == 0 ):

          dg = DELTAG_HI(rna,hs,js)

          dg = dg / 2.0e0

          x = beta * dg
          x = math.exp(-x) * rateh

          rna.wrk1[hs] = x
          rna.wrk1[js] = 0.0e0

          rna.wrk2[hs] = 0.0e0
          rna.wrk2[js] = 0.0e0

          atot += x

          hs += 1
          js -= 1

        #endwhile
      #endif

      if ( k != ke ):
        k = rna.ibsp[k]
        icnt += 1
      #endif

    #endif

    k += 1
    icnt += 1

  #endwhile

  rna.ptot[indx] = atot

  rna.LOOP_RESUM(indx)

  return
  