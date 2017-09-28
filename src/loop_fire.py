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

                                return

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

                    if ( iloop == 1 ) and ( nh == 2 and ns == 2 ):
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
                        rna.nsgl[indx] = ns -2

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

                        rna.nl -= 1

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

                    return

                #endif
            #endif

            #=== Helix Retraction ===#

            


            #=== Helix Morphing ===#




        #=== Defect Diffusion ===#




        #=== Open BP Inside Helix ===#