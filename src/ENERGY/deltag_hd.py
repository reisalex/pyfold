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

import math
from rnavar import em,eh,es,eaup,beta

def DELTAG_HD(rna,ii,jj,kk):

  # INTEGER
  # ii,jj,kk
  # i,j,k,l,n,ip,jp
  # nh,ns,mh,hs,js
  # indx,iloop,kloop
  #note replaced is with hs for Cython code

  # REAL
  # e1,e2,e3,e4,e5,ed
  # x,c,ei,ef

  indx = rna.link[jj]

  i = rna.loop[indx]
  j = rna.ibsp[i]
  n = rna.n

  nh = rna.nhlx[indx]
  nsh = rna.nsgl[indx]

  c = 1.750e0 / float(beta)

  #=== Internal Loop iloop = 1 ===#
  #=== External Loop iloop = 0 ===#

  iloop = 1
  kloop = 1

  if ( i > j ):
    iloop = 0
    i = 1
    j = n
  #endif

  if ( rna.link[ii] != 0 ):

    k = rna.link[ii]

    mh = rna.nhlx[k]
    ms = rna.nsgl[k]

    k = rna.loop[k]
    l = rna.ibsp[k]

    if ( k > l ):
      kloop = 0
      k = 1
      l = n
    #endif

  else:

    mh = 2
    ms = 0

    k = ii
    l = jj

  #endif

  e4 = 0.0e0

  #=== INITIAL ENERGY ===#

  if ( nh == 1 and iloop == 1 ):

    e1 = EHAIR(rna.iseq,i,j,n)

  elif ( nh == 2 and iloop == 1 ):

    ip = i + 1
    while ( rna.ibsp[ip] == 0 ):
      ip += 1
    #endwhile

    jp = rna.ibsp[ip]

    e1 = EBULGE(rna.iseq,i,j,ip,jp,n)

  else:

    e1 = 0.0e0
    e5 = 0.0e0
    e3 = 0.0e0

    hs = rna.iseq[ii]
    js = rna.iseq[jj]

    if ( iloop == 1 ):

      if ( ns <= 6 ):
        e1 = em + es * float(ns) + eh * float(nh)
      else:
        x = float(ns) / 6.0e0
        e1 = em + es * 6.0e0 + eh * float(nh)
        e1 += c * math.log(x)
      #endif
    #endif

    if ( ii > 1 ) and ( rna.ibsp[ii-1] == 0 ):

      e5 = EDANGLE(rna.iseq,ii,jj,ii-1,n)

      if ( ii > 2 ):

        ip = ii - 2
        jp = rna.ibsp[ip]

        if ( jp != 0 ):
          ed = EDANGLE(rna.iseq,ip,jp,ii-1,n)
          if ( kk == ii+1 ): e4 += ed
          e5 = min(e5,ed)
        #endif

      #endif

      if ( ii > 3 and kk == ii-1 ) and
         ( rna.ibsp[ii-2] == 0 ):

        ip = ii - 3
        jp = rna.ibsp[ip]

        if ( jp != 0 ):
          ed = EDANGLE(rna.iseq,ip,jp,ii-2,n)
          e1 += ed
        #endif
      #endif
    #endif

    if ( jj < n ) and ( rna.ibsp[jj+1] == 0 ):

      e3 = EDANGLE(rna.iseq,ii,jj,jj+1,n)

      if ( jj < n-1 ):

        ip = jj + 2
        jp = rna.ibsp[ip]

        if ( jp != 0 ):
          ed = EDANGLE(rna.iseq,ip,jp,jj+1,n)
          if ( kk == jj-1 ): e4 += ed
          e3 = min(e3,ed)
        #endif
      #endif

      if ( jj < n-2 and kk == jj+1 ) and
         ( rna.ibsp[jj+2] == 0 ):

        ip = jj + 3
        jp = rna.ibsp[ip]

        if ( jp != 0 ):
          ed = EDANGLE(rna.iseq,ip,jp,jj+2,n)
          e1 += ed
        #endif

      #endif
    #endif

    e1 = e1 + e3 + e5 + eaup[hs][js]

  #endif

  if ( mh == 1 and kloop == 1 ):

    e2 = EHAIR(rna.iseq,ii,jj,n)

  elif ( mh == 2 and kloop == 1 ):

    if ( ms == 0 ):

      e2 = ESTACK(rna.iseq,ii,jj,ii+1,jj-1,n)

    else:

      ip = ii + 1
      while ( rna.ibsp[ip] == 0 ):
        ip += 1
      #endwhile

      jp = rna.ibsp[ip]

      e2 = EBULGE(rna.iseq,ii,jj,ip,jp,n)

    #endif

  else:

    e2 = 0.0e0
    e5 = 0.0e0
    e3 = 0.0e0

    hs = rna.iseq[ii]
    js = rna.iseq[jj]

    if ( kloop == 1 ):

      if ( ms <= 6 ):
        e2 = em + es * float(ms) + eh * float(mh)
      else:
        x = float(ms) / 6.0e0
        e2 = em + es * 6.0e0 + eh * float(mh)
        e2 = e2 + c * math.log(x)
      #endif
    #endif

    if ( jj > 1 ) and ( rna.ibsp[jj-1] == 0 ):

      e5 = EDANGLE(rna.iseq,ii,jj,jj-1,n)

      if ( jj > 2 ):

        ip = jj - 2
        jp = rna.ibsp[ip]

        if ( jp != 0 ):
          ed = EDANGLE(rna.iseq,ip,jp,jj-1,n)
          if ( kk == jj+1 ): e4 += ed
          e5 = min(e5,ed)
        #endif
      #endif

      if ( jj > 3 and kk == jj-1 ) and ( rna.iseq[jj-2] == 0 ):

        ip = jj - 3
        jp = rna.ibsp[ip]

        if ( jp != 0 ):
          ed = EDANGLE(rna.iseq,ip,jp,jj-2,n)
          e2 += ed
        #endif
      #endif
    #endif

    if ( ii < n ) and ( rna.ibsp[ii+1] == 0 ):

      e3 = EDANGLE(rna.iseq,ii,jj,ii+1,n)

      if ( ii < n-1 ):

        ip = ii + 2
        jp = rna.ibsp[ip]

        if ( jp != 0 ):
          ed = EDANGLE(rna.iseq,ip,jp,ii+1,n)
          if ( kk == ii-1 ): e4 += ed
          e3 = min(e3,ed)
        #endif
      #endif

      if ( ii < n-2 and kk == ii+1 ) and
         ( rna.ibsp[ii+2] == 0 ):

        ip = ii + 3
        jp = rna.ibsp[ip]

        if ( jp != 0 ):
          ed = EDANGLE(rna.iseq,ip,jp,ii+2,n)
          e2 += ed
        #endif
      #endif
    #endif

    e2 = e2 + e3 + e5 + eaup[hs][js]

  #endif

  ei = e1 + e2

  #=== FINAL ENERGY ===#

  if ( kk == ii-1 or kk == jj+1 ):
    ns -= 1
    ms += 1
  else:
    ns += 1
    ms -= 1
  #endif

  if ( nh == 1 and iloop == 1 ):

    if ( kk == i-1 or kk == i+1 ):
      e1 = EHAIR(rna.iseq,kk,j,n)
    #endif

    if ( kk == j-1 or kk == j+1 ):
      e1 = EHAIR(rna.iseq,i,kk,n)
    #endif

  elif ( nh == 2 and iloop == 1 ):

    if ( ns == 0 ):

      if ( kk == ii-1 ):
        e1 = ESTACK(rna.iseq,kk-1,jj+1,kk,jj,n)
      else:
        e1 = ESTACK(rna.iseq,ii-1,kk+1,ii,kk,n)
      #endif

    else:

      ip = i + 1
      while ( rna.ibsp[ip] == 0 ):
        ip += 1
      #endwhile

      jp = rna.ibsp[ip]

      if ( ii == j ):

        if ( kk == ii+1 or kk == ii-1 ):
          e1 = EBULGE(rna.iseq,jj,kk,ip,jp,n)
        else:
          e1 = EBULGE(rna.iseq,kk,ii,ip,jp,n)
        #endif

      else:

        if ( kk == ii+1 or kk == ii-1 ):
          e1 = EBULGE(rna.iseq,i,j,kk,jj,n)
        else:
          e1 = EBULGE(rna.iseq,i,j,ii,kk,n)
        #endif
      #endif

    #endif

  else:

    e1 = 0.0e0
    e5 = 0.0e0
    e3 = 0.0e0

    if ( iloop == 1 ):

      if ( ns <= 6 ):
        e1 = em + es * float(ns) + eh * float(nh)
      else:
        x = float(ns) / 6.0e0
        e1 = em + es * 6.0e0 + eh * float(nh)
        e1 = e1 + c * math.log(x)
      #endif
    #endif

    if ( kk == jj+1 or kk == jj-1 ):

      hs = rna.iseq[ii]
      js = rna.iseq[kk]

      if ( ii > 1 ) and ( rna.ibsp[ii-1] == 0 ):

        e5 = EDANGLE(rna.iseq,ii,kk,ii-1,n)

        if ( ii > 2 ):

          ip = ii - 2
          jp = rna.ibsp[ip]

          if ( jp != 0 ):
            ed = EDANGLE(rna.iseq,ip,jp,ii-1,n)
            e5 = min(e5,ed)
          #endif
        #endif
      #endif

      if ( kk < n ) and ( rna.ibsp[kk+1] == 0 or kk == jj-1 ):

        e3 = EDANGLE(rna.iseq,ii,kk,kk+1,n)

        if ( kk < n-1 ):

          ip = kk + 2
          ip = rna.ibsp[ip]

          if ( jp != 0 ):
            ed = EDANGLE(rna.iseq,ip,jp,kk+1,n)
            e3 = min(e3,ed)
          #endif
        #endif
      #endif

    #endif
    
    if ( kk == ii-1 or kk == ii+1 ):

      hs = rna.iseq[kk]
      js = rna.iseq[jj]

      if ( kk > 1 ) and (rna.ibsp[kk-1] == 0 or kk == ii+1 ):

        e5 = EDANGLE(rna.iseq,kk,jj,kk-1,n)

        if ( kk > 2 ):

          ip = kk - 2
          jp = rna.ibsp[ip]

          if ( jp != 0 ):
            ed = EDANGLE(rna.iseq,ip,jp,kk-1,n)
            e5 = min(e5,ed)
          #endif
        #endif

      #endif

      if ( jj < n ) and ( rna.ibsp[jj+1] == 0 ):

        e3 = EDANGLE(rna.iseq,kk,jj,jj+1,n)

        if ( jj < n - 1 ):

          ip = jj + 2
          jp = rna.ibsp[ip]

          if ( jp != 0 ):
            ed = EDANGLE(rna.iseq,ip,jp,jj+1,n,ed)
            e3 = min(e3,ed)
          #endif
        #endif
      #endif

    #endif

    e1 = e1 + e3 + e5 + eaup[hs][js]

  #endif

  if ( mh == 1 and kloop == 1 ):

    if ( kk == ii-1 or kk == ii+1 ):
      e2 = EHAIR(rna.iseq,kk,jj,n)
    #endif

    if ( kk == jj-1 or kk == jj+1 ):
      e2 = EHAIR(rna.iseq,ii,kk,n)
    #endif

  elif ( mh == 2 and kloop == 1 ):

    if ( ms == 0 ):

      if ( kk == ii+1 ):
        e2 = ESTACK(rna.iseq,kk,jj,kk+1,jj-1,n)
      else:
        e2 = ESTACK(rna.iseq,ii,kk,ii+1,kk-1,n)
      #endif
    
    else:

      ip = ii + 1
      while ( rna.ibsp[ip] == 0 ):
        ip += 1
      #endwhile

      jp = rna.ibsp[ip]

      if ( kk == ii+1 or kk == ii-1 ):
        e2 = EBUGLE(rna.iseq,kk,jj,ip,jp,n)
      else:
        e2 = EBUGLE(rna.iseq,ii,kk,ip,jp,n)
      #endif
    #endif

  else:

    e2 = 0.0e0
    e5 = 0.0e0
    e3 = 0.0e0

    if ( kloop == 1 ):

      if ( ms <= 6 ):
        e2 = em + es * float(ms) + eh * float(mh)
      else:
        x = float(ms) / 6.0e0
        e2 = em + es * 6.0e0 + eh * float(mh)
        e2 = e2 + c * math.log(x)
      #ndif
    #endif

    if ( kk == jj+1 or kk == jj-1 ):

      hs = rna.iseq[ii]
      js = rna.iseq[kk]

      if ( kk > 1 ) and ( rna.ibsp[kk-1] == 0 or kk == jj+1 ):

        e5 = EDANGLE(rna.iseq,ii,kk,kk-1,n)

        if ( kk > 2 ):

          ip = kk - 2
          jp = rna.ibsp[ip]

          if ( jp != 0 ):
            ed = EDANGLE(rna.iseq,ip,jp,kk-1,n)
            e5 = min(e5,ed)
          #endif
        #endif
      #endif

      if ( ii < n ) and ( rna.ibsp[ii+1] == 0 ):

        e3 = EDANGLE(rna.iseq,ii,kk,ii+1,n)

        if ( ii < n-1 ):

          ip = ii + 2
          jp = rna.ibsp[ip]

          if ( jp != 0 ):
            ed = EDANGLE(rna.iseq,ip,jp,ii+1,n)
            e3 = min(e3,ed)
          #endif
        #endif
      #endif

    #endif

    if ( kk == ii-1 or kk == ii+1 ):

      hs = rna.iseq[kk]
      js = rna.iseq[jj]

      if ( jj > 1 ) and ( rna.ibsp[jj-1] == 0 ):

        e5 = EDANGLE(rna.iseq,kk,jj,jj-1,n)

        if ( jj > 2 ):

          ip = jj - 2
          jp = rna.ibsp[ip]

          if ( jp != 0 ):
            ed = EDANGLE(rna.iseq,ip,jp,jj-1,n)
            e5 = min(e5,ed)
          #endif
        #endif
      #endif

      if ( kk < n ) and ( rna.ibsp[kk+1] == 0 or kk == ii-1 ):

        e3 = EDANGLE(rna.iseq,kk,jj,kk+1,n)

        if ( kk < n-1 ):

          ip = kk + 2
          jp = rna.ibsp[ip]

          if ( jp != 0 ):
            ed = EDANGLE(rna.iseq,ip,jp,kk+1,n)
            e3 = min(e3,ed)
          #endif
        #endif
      #endif
    #endif

    e2 = e2 + e3 + e5 + eaup[hs][js]

  #endif

  ef = e1 + e2 + e4

  dg = float(ef) - float(ei)

  return dg





