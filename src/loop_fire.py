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

def LOOP_FIRE(rna,indx,amax):

  # FLOAT
  # amax,x,dg,atot

  # INTEGER
  # idnx
  # i,j,k,n,nl,ip,jp,kp,hs,js
  # ks,ke,l,lmx,icnt,iloop
  # nt,nh,ns,mt,mh,ms,icase
  # jndx,kndx,nsum

  # the FORTRAN <is> variable was converted to <hs> for the Python code

  i = rna.loop[indx]
  j = rna.ibsp[i]

  n = rna.n
  nl = rna.nl

  nsum = rna.nsum

  if ( i == n ): j = 1

  nh = rna.nhlx[indx]
  ns = rna.nsgl[indx]

  #=== Internal Loop iloop = 1 ===#
  #=== External Loop iloop = 0 ===#

  if ( i < j ): iloop = 1
  if ( i > j ): iloop = 0

  if ( iloop == 1 ):
    ks = i + 1
    ke = j
  else:
    ks = j
    ke = i

  #=== FIND REACTION TO FIRE ===#
  
  atot = 0.0e0
  jndx = 0
  kndx = 0

  k = ks
  icnt = 0

  while ( k <= ke ):

    #=== Nucleation Events ===#

    if ( rna.ibsp[k] == 0 ):

      if ( atot + rna.wrk1[k] >= amax ):

        l = 2
        lmx = nt / 2 + 1

        if ( nt % 2 == 0 ) and ( icnt+1 > lmx-1 ): lmx -= 1

        if ( iloop == 0 ): lmx = nt - icnt

        kp = k + 1
        hs = rna.iseq[k]

        while ( l <= lmx ):

          if ( rna.ibsp[kp] == 0 ):

            js = rna.iseq[kp]

            if ( l > 4 and iwc[hs][js] == 1 ):

              atot += pnuc[l]

              if ( atot >= amax ):

                nl += 1

                rna.ibsp[k] = kp
                rna.ibsp[kp] = k
                rna.nl = nl

                if ( nl > nsum ):
                  rna.nsum = 2 * nsum

                if ( k < kp ):
                  rna.loop[nl] = k
                  rna.link[k]  = nl
                  rna.link[kp] = indx
                else:
                  rna.loop[nl] = kp
                  rna.link[k]  = indx
                  rna.link[kp] = nl

                # Fix links in new loop

                mh = 1
                ms = 0

                ip = min(k,kp)
                jp = ip + 1

                jndx = rna.link[ip]

                while ( jp < rna.ibsp[ip] ):

                  if ( rna.link[jp] == indx ):
                    rna.link[jp] = jndx

                  if ( rna.ibsp[jp] > jp ): mh += 1
                  if ( rna.ibsp[jp] == 0 ): ms += 1

                  if ( rna.ibsp[jp] > jp ):
                    jp = rna.ibsp[jp]
                  else:
                    jp += 1

                rna.nhlx[indx] = nh - mh + 2
                rna.nsgl[indx] = ns - ms - 2

                rna.nhlx[jndx] = mh
                rna.nsgl[jndx] = ms

                rna.LOOP_REAC(indx)
                rna.LOOP_REAC(jndx)

                # Recalc lower loop?

                if ( iloop == 0 ): kndx = 0
                if ( iloop == 1 ): kndx = rna.link[j]

                if ( kndx != 0 ): rna.LOOP_REAC(kndx)

                return rna

          else:

            l += 1
            kp = rna.ibsp[kp]

          l += 1
          kp += 1

      atot += rna.wrk1[k]

    #=== Helix Events ===#

    if ( rna.ibsp[k] > 0 ):

      ip = k
      jp = rna.ibsp[k]

      icase = 0

      #=== Helix Extension ===#

      if ( ip > 1 and jp < n ) and ( nh > 1 or ns > 4 ) and \
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

      if ( icase > 0 ):

        atot += rna.wrk2[ip]

        if ( atot >= amax ):

          if ( icase == 1 ):

            jndx = rna.link[ip]
            rna.nsgl[indx] = ns - 2

          elif ( icase == 2 ):

            jndx = rna.link[i+2]

            if ( jndx == nl ): jndx = indx

            # Delete loop indx
            # Copy loop nl to indx

            if ( indx != nl ):

              rna.loop[indx] = rna.loop[nl]
              rna.nhlx[indx] = rna.nhlx[nl]
              rna.nsgl[indx] = rna.nsgl[nl]
              rna.ptot[indx] = rna.ptot[nl]

              rna.LOOP_RESUM(indx)

              kp = rna.loop(nl)

              rna.link[kp] = indx

              l = kp + 1

              while ( l < rna.ibsp[kp] ):

                if ( rna.link[l] == nl ):
                  rna.link[l] = indx

                if ( rna.ibsp[l] > l ):
                  l = rna.ibsp[l]
                else:
                  l += 1
                # endif
              #endwhile
            #endif

            rna.link[i]   = 0
            rna.link[j-2] = 0

            rna.loop[nl]  = 0
            rna.nhlx[nl]  = 0
            rna.nsgl[nl]  = 0
            rna.ptot[nl]  = 0.0e0

            rna.LOOP_RESUM(nl)

            rna.nl = nl - 1

            if ( nsum > 2 ) and ( rna.nl <= nsum / 2 ):
              nsum = nsum / 2
              rna.nsum = nsum
              rna.psum[nsum] = 0.0e0
            #endif

          elif ( icase == 3 ):

            rna.loop[indx] = jp + 1
            rna.nsgl[indx] = ns - 2

          #endif

          # Adjust base-pairs
          rna.ibsp[ip-1] = jp + 1
          rna.ibsp[jp+1] = ip - 1

          # Fix links
          rna.link[jp+1] = rna.link[jp]
          rna.link[jp] = 0

          # Recalc main loop?
          if ( icase != 2 ):
            rna.LOOP_REAC(indx)
          #endif

          # Recalc upper loop?
          if ( icase != 3 ):
            if ( jndx == 0 ):
              rna.HELX_REAC(ip)
            else:
              rna.LOOP_REAC(jndx)
            #endif
          #endif

          # Recalc lower loop?
          if ( iloop == 0 ): kndx = 0
          if ( iloop == 1 ): kndx = rna.link[j]

          if ( kndx != 0 ): rna.LOOP_REAC(kndx)

          return rna

        #endif
      #endif

      #=== Helix Retraction ===#

      icase = 0

      if ( ip != n and jp != 1 ):

        if ( rna.ibsp[ip+1] == jp-1 ):

          icase = 1

          if ( iloop == 1 ) and ( k == ke ):
            icase = 2
          #endif

        elif ( iloop == 0 or k != ke ):

          l = rna.link[ip]
          mh = rna.nhlx[l]
          ms = rna.nsgl[l]

          icase = 3

        #endif
      #endif

      if ( icase > 0 ):

        if ( icase == 2 ):
          atot += rna.wrk1[ip]
        else:
          atot += rna.wrk1[jp]
        #endif

        if ( atot >= amax ):

          rna.ibsp[ip] = 0
          rna.ibsp[jp] = 0

          if ( icase == 1 ):

            rna.nsgl[indx] = ns + 2

            rna.link[jp-1] = rna.link[jp]
            if ( jp != n ): rna.link[jp] = 0

            # Recalc main loop

            rna.LOOP_REAC(indx)

            # Recalc upper/lower loops?

            jndx = rna.link[ip+1]

            if ( iloop == 0 ): kndx = 0
            if ( iloop == 1 ): kndx = rna.link[j]

            if ( jndx == 0 ): rna.HELX_REAC(ip+1)
            if ( jndx != 0 ): rna.LOOP_REAC(jndx)
            if ( kndx != 0 ): rna.LOOP_REAC(kndx)

          elif ( icase == 2 ):

            rna.nsgl[indx] = ns + 2
            rna.loop[indx] = jp - 1

            rna.link[jp-1] = rna.link[jp]
            if ( jp != n ): rna.link[jp] = 0

            # Recalc main loop

            rna.LOOP_REAC(indx)

            # Recalc lower loop?

            kndx = rna.link[ip+1]

            if ( kndx != 0 ): rna.LOOP_REAC(kndx)

          elif ( icase == 3 ):

            rna.nhlx[indx] = nh + mh - 2
            rna.nsgl[indx] = ns + ms + 2

            jndx = rna.link[ip]

            # Fix links in loop being deleted

            l = ip + 1

            while ( l < jp ):

              if ( rna.link[l] == jndx ):
                rna.link[l] = indx

              if ( rna.ibsp[l] > l ):
                l = rna.ibsp[l]
              else:
                l += 1
              #endif
            #endwhile

            # Copy loop nl to jndx

            if ( jndx != nl ):

              rna.loop[jndx] = rna.loop[nl]
              rna.nhlx[jndx] = rna.nhlx[nl]
              rna.nsgl[jndx] = rna.nsgl[nl]
              rna.ptot[jndx] = rna.ptot[nl]

              rna.LOOP_RESUM(jndx)

              kp = rna.loop[nl]

              rna.link[kp] = jndx

              l = kp + 1

              while ( l < rna.ibsp[kp] ):

                if ( rna.link[l] == nl ):
                  rna.link[l] = jndx
                #endif

                if ( rna.ibsp[l] > l ):
                  l = rna.ibsp[l]
                else:
                  l += 1
                #endif

              #endwhile
            #endif

            if ( indx != nl ): jndx = indx

            rna.link[ip] = 0
            if ( jp != n ): rna.link[jp] = 0

            rna.loop[nl] = 0
            rna.nhlx[nl] = 0
            rna.nsgl[nl] = 0
            rna.ptot[nl] = 0.0e0

            rna.LOOP_RESUM(nl)

            rna.nl = nl - 1

            if ( nsum > 2 ) and ( rna.nl <= nsum / 2 ):
              nsum = nsum / 2
              rna.nsum = nsum
              rna.psum[nsum] = 0.0e0
            #endif

            # Recalc main loop

            rna.LOOP_REAC(jndx)

            # Recalc lower loop?

            if ( iloop == 0 ): kndx = 0
            if ( iloop == 1 ): kndx = rna.link[j]

            if ( kndx != 0 ): rna.LOOP_REAC(kndx)

          #endif

          return rna

        #endif
      #endif

      icase = 0

      #=== Helix Morphing ===#

      if ( iloop == 0 or nh > 2 ) and (ip > 1 and jp < n):

        hs = rna.iseq[ip-1]
        js = rna.iseq[jp+1]

        if ( iwc[hs][js] == 1 ): icase = 1

        hs = rna.ibsp[ip-1]
        js = rna.ibsp[jp+1]

        if ( hs != 0 ) and ( rna.link[hs] != 0 ): icase = 0

        if (js != 0 ) and (rna.link[jp+1] != 0 ): icase = 0

        if ( hs == 0 and js == 0 ): icase = 0

      #endif

      if ( icase > 0 ):

        dg = DELTAG_HM(rna,ip,jp)

        dg = dg / 2.0e0

        x = beta * dg
        x = math.exp(-x) * ratem

        atot += x

        if ( atot >= amax ):

          hs = rna.ibsp[ip-1]
          js = rna.ibsp[jp+1]

          if ( hs != 0 ):

            rna.ibsp[hs]   = 0
            rna.link[ip-1] = 0
            rna.link[ip-2] = indx

            if ( iloop == 1 and hs == ke ):
              rna.loop[indx] = ip - 2
            #endif
          #endif

          if ( js != 0 ):

            rna.ibsp[js] = 0
            rna.link[js] = 0
            rna.link[js-1] = indx

            if ( iloop == 1 and js == i ):
              rna.loop[indx] = js - 1
            #endif
          #endif

          if ( hs != 0 and js != 0 ):
            
            ns += 2
            rna.nsgl[indx] = ns
          #endif

          # Adjust base pairs

          rna.ibsp[ip-1] = jp+1
          rna.ibsp[jp+1] = ip-1

          # Fix links

          rna.link[jp+1] = rna.link[jp]
          rna.link[jp] = 0

          if ( iloop == 1 and k = ke ):
            rna.loop[indx] = jp + 1
          #endif

          # Recalc main loop

          rna.LOOP_REAC(indx)

          jndx = rna.link[ip]
          if ( jndx == 0 ): rna.HELX_REAC(ip)
          else:             rna.LOOP_REAC(jndx)

          # Recalc 5' 3' loops?

          if ( hs != 0 ):
            jndx = rna.link[hs+1]
            if ( jndx == 0 ): rna.HELX_REAC(hs+1)
            else:             rna.LOOP_REAC(jndx)
          #endif

          if ( js != 0 ):
            jndx = rna.link[jp+2]
            if ( jndx == 0 ): rna.HELX_REAC(jp+2)
            else:             rna.LOOP_REAC(jndx)
          #endif

          # Recalc lower loop?

          if ( iloop == 0 ): kndx = 0
          if ( iloop == 1 ): kndx = rna.link[j]

          if ( kndx != 0 ): rna.LOOP_REAC(kndx)

          return rna

        #endif
      #endif

      #=== Defect Diffusion ===#

      # PUSH

      icase = 0

      if ( rna.link[ip] == 0 ):

        icase = 2

        if ( iloop == 1 ):
          if ( nh == 2 and ns == 1 ): icase = 3
          if ( nh == 1 and ns == 3 ): icase = 0
        #endif

      elif ( iloop == 0 or k != ke ):

        icase = 1

        if ( iloop == 1 ) and ( nh == 2 and ns == 1):
          icase = 4
        #endif

      #endif

      if ( icase > 0 ):

        # PUSH 5' end

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

            if ( atot >= amax ):
              return GOTO(rna)

          #endif
        #endif

        # PUSH 3' end

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

            if ( atot > amax ):
              return GOTO(rna)

          #endif
        #endif
      #endif

      # PULL

      icase = 0

      if ( rna.link[ip] != 0 ) and ( iloop == 0 or k != ke ):

        l  = rna.link[ip]
        mh = rna.nhlx[l]
        ms = rna.nsgl[l]

        icase = 1

        if ( mh == 1 and ms == 3 ): icase = 0
        if ( mh == 2 and ms == 1 ): icase = 4

      #endif

      if ( icase > 0 ):

        # PULL 5' end

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

            if ( atot >= amax ):
              return GOTO(rna)

          #endif
        #endif

        # PULL 3' end

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

            if ( atot >= amax ):
              return GOTO(rna)

          #endif
        #endif
      #endif

      if ( atot >= amax ):
        return GOTO(rna)
      #endif

      #=== Open BP Inside Helix ===#

      if ( rna.link[ip] == 0 ) and ( iloop == 1 and k == ke ):

        hs = ip + 1
        js = jp - 1

        while ( rna.link[hs] == 0 ):

          atot += rna.wrk1[hs]

          hs += 1
          js -= 1

          if ( atot >= amax ):

            nl += 1

            rna.ibsp[hs-1] = 0
            rna.ibsp[js+1] = 0
            rna.nl = nl

            rna.loop[nl] = js
            rna.link[js] = nl
            rna.link[hs-2] = nl

            rna.nhlx[nl] = 2
            rna.nsgl[nl] = 2

            if ( nl > nsum ):
              rna.nsum = 2 * nsum
            #endif

            # Recalc upper loop?

            jndx = rna.link[js+2]

            if ( jndx == 0 ): rna.HELX_REAC(js+2)
            else:             rna.LOOP_REAC(jndx)

            # Calculate new loop

            jndx = rna.link[js]

            rna.LOOP_REAC(jndx)

            # Recalc lower loop?

            kndx = rna.link[hs]

            if ( kndx != 0 ): rna.LOOP_REAC(kndx)

            return rna

          #endif
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

  return rna


# Replacing the GOTO structure of KFOLD's loop_fire.f90 subroutine

def GOTO(rna): # Going to need more args

  # Adjust base pairs

  if ( kp == ip+1 or kp == ip-1 ):
    rna.ibsp[ip] = 0
    rna.ibsp[jp] = kp
    rna.ibsp[kp] = jp
  else:
    rna.ibsp[ip] = kp
    rna.ibsp[jp] = 0
    rna.ibsp[kp] = ip
  #endif

  if ( icase == 1 ):

    jndx = rna.link[ip]

    if ( kp == ip+1 or kp == ip-1 ):
      rna.loop[jndx] = kp
      rna.link[kp] = rna.link[ip]
      rna.link[ip] = 0
    else:
      rna.link[kp] = rna.link[jp]
      rna.link[jp] = 0
    #endif

    if ( kp == ip+1 or kp == jp-1 ):
      rna.nsgl[indx] += 1
      rna.nsgl[jndx] -= 1
    else:
      rna.nsgl[indx] -= 1
      rna.nsgl[jndx] += 1
    #endif

    rna.LOOP_REAC(indx)
    rna.LOOP_REAC(jndx)

    # Recalc lower loop?

    if ( iloop == 0 ): kndx = 0
    if ( iloop == 1 ): kndx = rna.link[j]

    if ( kndx != 0 ): rna.LOOP_REAC(kndx)

  elif ( icase == 2 ):

    # Add loop

    nl += 1

    rna.nl = nl

    r.link[jp-1] = nl

    if ( nl > nsum ):
      rna.nsum = 2 * nsum
    #endif

    if ( kp == jp+1 ):
      rna.link[ip] = nl
      rna.link[jp+1] = rna.link[jp]
      rna.link[jp] = 0
    else:
      rna.link[kp] = nl
    #endif

    if ( iloop == 1 and k == ke ):
      rna.loop[nl] = i - 1
      if ( kp == jp+1 ):
        rna.loop[indx] = i + 1
      #endif
    elif ( kp == ip-1 ):
      rna.loop[nl] = kp
    else:
      rna.loop[nl] = ip
    #endif

    rna.nsgl[indx] = ns - 1
    rna.nhlx[nl] = 2
    rna.nsgl[nl] = 1

    rna.LOOP_REAC(indx)
    rna.LOOP_REAC(nl)

    # Recalc upper loop?

    jndx = rna.link[ip+1]
    if ( jndx != 0 ): rna.LOOP_REAC(jndx)

    if ( iloop == 0 or k != ke ) and (jndx == 0 ):
      rna.HELX_REAC(ip+1)
    #endif

    # Recalc lower loop?

    if ( iloop == 0 ): kndx = 0
    if ( iloop == 1 ): kndx = rna.link[j]

    if ( kndx != 0 ): rna.LOOP_REAC(kndx)

  elif ( icase == 3 ):

    # Move loop

    rna.link[jp-1] = rna.link[jp]
    rna.link[jp] = 0

    if ( kp == jp+1 ):
      rna.link[ip] = rna.link[ip-1]
      rna.link[ip-1] = 0
    else:
      rna.link[ip-1] = rna.link[ip-2]
      rna.link[ip-2] = 0
    #endif

    if ( k == ke ):

      rna.loop[indx] = i-1

      kndx = rna.link[ip+1]

      if ( kp == ip-1 ): jndx = rna.link[i+1]
      if ( kp == jp+1 ): jndx = rna.link[i+2]

      if ( jndx == 0 ): rna.HELX_REAC(j-1)

    else:

      rna.loop[indx] = i+1

      jndx = rna.link[ip+1]
      kndx = rna.link[j]

      if ( jndx == 0 ): rna.HELX_REAC(ip+1)

    #endif

    rna.LOOP_REAC(indx)

    if ( jndx != 0 ): rna.LOOP_REAC(jndx)
    if ( kndx != 0 ): rna.LOOP_REAC(kndx)

  elif ( icase == 4 ):

    # Delete loop

    if ( kp == ip+1 ):

      jndx = rna.link[ip]

      rna.nsgl[indx] += 1

      rna.link[ip] = 0
      rna.link[jp-1] = 0

      rna.LOOP_REAC(indx)

      # Recalc upper loop?

      kndx = rna.link[ip+2]

      if ( kndx != 0 ): rna.LOOP_REAC(kndx)
      if ( kndx == 0 ): rna.HELX_REAC(r,ip+2)

      # Recalc lower loop?

      if ( iloop == 0 ): kndx = 0
      if ( iloop == 1 ): kndx = rna.link[j]

      if ( kndx != 0 ): rna.LOOP_REAC(kndx)

    elif ( kp == jp-1 ):

      jndx = rna.link[ip]

      rna.nsgl[indx] += 1

      rna.link[jp-1] = rna.link[jp]
      rna.link[ip]   = 0
      rna.link[jp]   = 0
      rna.link[jp-2] = 0

      rna.LOOP_REAC(indx)

      # Recalc upper loop?

      kndx = rna.link[ip+1]

      if ( kndx != 0 ): rna.LOOP_REAC(kndx)
      if ( kndx == 0 ): rna.HELX_REAC(ip+1)

      # Recalc lower loop?

      if ( iloop == 0 ): kndx = 0
      if ( iloop == 1 ): kndx = rna.link[j]

      if ( kndx != 0 ): rna.LOOP_REAC(kndx)

    else:

      jndx = rna.link[ip]
      kndx = rna.link[j]

      rna.nsgl[jndx] += 1

      rna.link[i]  = 0
      rna.link[jp] = 0

      if ( kp == ip-1 ):
        rna.loop[jndx] = ip - 1
        rna.link[ip-1] = rna.link[ip]
        rna.link[ip] = 0
      #endif

      # Recalc loops

      rna.LOOP_REAC(jndx)

      if ( kndx != 0 ): rna.LOOP_REAC(kndx)

      jndx = indx

    #endif

    # Delete loop jndx

    # Copy loop nl to jndx

    if ( jndx != nl ):

      rna.loop[jndx] = rna.loop[nl]
      rna.nhlx[jndx] = rna.nhlx[nl]
      rna.nsgl[jndx] = rna.nsgl[nl]
      rna.ptot[jndx] = rna.ptot[nl]

      rna.LOOP_RESUM(jndx)

      hs = rna.loop[nl]

      rna.link[hs] = jndx

      l = hs + 1

      while ( l < rna.ibsp[hs] ):

        if ( rna.link[l] == nl ):
          rna.link[l] = jndx
        #endif

        if ( rna.ibsp[l] > l ):
          l = rna.ibsp[l]
        else:
          l += 1
        #endif

      #endwhile

    #endif

    rna.loop[nl] = 0
    rna.nhlx[nl] = 0
    rna.nsgl[nl] = 0
    rna.ptot[nl] = 0.0e0

    rna.LOOP_RESUM(nl)

    rna.nl = nl - 1

    if ( nsum > 2 ) and ( rna.nl <= nsum / 2 ):
      nsum = nsum / 2
      rna.nsum = nsum
      rna.psum[nsum] = 0.0e0
    #endif
  #endif

  return rna