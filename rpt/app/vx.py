#!/usr/bin/env python

'''
corrections for the three estimators
0th order for < exp(-beta e) > = 1 
1st order for < e min{ exp(-beta e), 1} > = 0
2nd order for < e exp(-beta/2 e) > = 0
'''

import os, sys
from math import *

#fnprof = "profis16.dat" # profile name
#fndos = "../islogdos16x16.txt" # density of states file
fnprof = "profis32.dat" # profile name
fndos = "../islogdos32x32.txt" # density of states file
s = open(fnprof, "r").readlines()


n = len(s)
ene = [0]*n
betx = [0]*n
bets = [0]*n
beth = [0]*n
refbet = [0]*n
serrx = [0]*n
serrs = [0]*n
serrh = [0]*n
hn8 =[0]*n
hn4 =[0]*n
h0 =[0]*n
hp4 =[0]*n
hp8 =[0]*n
crx1 = [0]*n
crx2 = [0]*n
crs1 = [0]*n
crs2 = [0]*n
crh1 = [0]*n
crh2 = [0]*n
mmx1 = [0]*n
mmx2 = [0]*n
mms2 = [0]*n
mms3 = [0]*n
mmh2 = [0]*n
lndos = [0]*n


# integrate bet 
sbetx = [0]*n
sbets = [0]*n
sbeth = [0]*n

s2 = open(fndos, "r").readlines()
for i in range(n):
  lndos[i] = float(s2[i].strip()[:16])

for i in range(n):
  arr = s[i].strip().split()
  ene[i]    = float(arr[0])
  betx[i]   = float(arr[4])
  bets[i]   = float(arr[6])
  beth[i]   = float(arr[8])
  refbet[i] = float(arr[9])
  serrx[i]   = float(arr[10])
  serrs[i]   = float(arr[12])
  serrh[i]   = float(arr[14])
  hn8[i] = float(arr[15])
  hn4[i] = float(arr[16])
  h0[i]  = float(arr[17])
  hp4[i] = float(arr[18])
  hp8[i] = float(arr[19])
  tot = hn8[i] + hn4[i] + h0[i] + hp4[i] + hp8[i]
  
  if tot > 0.:
    mmx2[i] = (hn8[i]*64.0 + hn4[i]*16.0 + hp4[i]*16.0 + hp8[i]*64.0)/tot
    mmx1[i] = (-hn8[i]*8.0 - hn4[i]*4.0 + hp4[i]*4.0 + hp8[i]*8.0)/tot

    if bets[i] > 0:
      xp4b = exp(-bets[i]*4)
      xp8b = exp(-bets[i]*8)
      xpn4b = 1
      xpn8b = 1
    else:
      xpn4b = exp(bets[i]*4)
      xpn8b = exp(bets[i]*8)
      xp4b = 1
      xp8b = 1

    mms2[i] = (hn8[i]*xpn8b*8.0**2 + hn4[i]*xpn4b*4.0**2 + hp4[i]*xp4b*4.0**2 + hp8[i]*xp8b*8.0**2)/tot
    mms3[i] = (hn8[i]*xpn8b*8.0**3 + hn4[i]*xpn4b*4.0**3 + hp4[i]*xp4b*4.0**3 + hp8[i]*xp8b*8.0**3)/tot
    if bets[i] > 0: mms3[i] *= -1
    
    '''
    if i < n/2:
      mms2[i] = (hn8[i]*64.0 + hn4[i]*16.0)/tot
      mms3[i] = -(hn8[i]*8.0**3 + hn4[i]*4.0**3)/tot
    elif i > n/2:
      mms2[i] = (hp8[i]*64.0 + hp4[i]*16.0)/tot
      mms3[i] = (hp8[i]*8.0**3 + hp4[i]*4.0**3)/tot
    '''
    
    
    if beth[i] > 100: beth[i] = 100
    if beth[i] < -100: beth[i] = -100
    xph4b = exp(0.5*beth[i]*4)
    xph8b = exp(0.5*beth[i]*8)

    mmh2[i] = (hn8[i]*xph8b*64.0 + hn4[i]*xph4b*16.0 + hp4[i]*1/xph4b*16.0 + hp8[i]*1/xph8b*64.0 )/tot


#
# I. Corrections for  < exp (- beta e) > = 1 
#

# recompute the density of states
sbetx[0] = log(2)
for i in range(1, n):
  ds = .5*4*(betx[i-1] + betx[i])  # smooth it
  if abs(ds) > 100: ds = 0
  sbetx[i] = sbetx[i-1] + ds

# compute the energy differential
dbdux = [0]*n
for i in range(n):
  if i == 0 or i == n - 1:
    continue
  dbdux[i] = (betx[i+1] - betx[i-1])/8

# the first correction term
# \int (1/2) dbdu < e^2 > / < e > 
# integrate from the bottom
for i in range(1, n/2+1):
  if mmx1[i-1] == 0 or mmx1[i] == 0:
    dcr1 = 0
  else:
    dcr1 = (dbdux[i-1]*mmx2[i-1] + dbdux[i]*mmx2[i])/(mmx1[i-1] + mmx1[i])
  crx1[i] = crx1[i-1] + .5*dcr1*4
cr1mida = crx1[n/2]

# integrate from the top
for i in range(n-1, n/2, -1):
  if mmx1[i-1] == 0 or mmx1[i] == 0:
    dcr1 = 0
  else:
    dcr1 = (dbdux[i-1]*mmx2[i-1] + dbdux[i]*mmx2[i])/(mmx1[i-1] + mmx1[i])
  crx1[i-1] = crx1[i] - .5*dcr1*4

cr1midb = crx1[n/2]
crx1[n/2] = .5 * (cr1mida + cr1midb)

# the second correction term
# - \int (log)' < e >  = - log < e >
for i in range(1, n):
  if fabs(mmx1[i]) > 0:
    if i < n/2:
      crx2[i] = -log(mmx1[i]/mmx1[n/2-1])
    elif i > n/2:
      crx2[i] = -log(mmx1[i]/mmx1[n/2+1])

#
# Corrections for < e min {exp(-bet * e), 1} > = 0
#

# recompute the density of states
sbets[0] = log(2)
for i in range(1, n):
  ds = .5*4*(bets[i-1] + bets[i])
  if abs(ds) > 100: ds = 0
  sbets[i] = sbets[i-1] + ds

# compute the energy differential
dbdus = [0]*n
for i in range(n):
  if i == 0 or i == n - 1:
    continue
  dbdus[i] = (bets[i+1] - bets[i-1])/8

# the first correction term
# \int (1/2) dbdu < e^3 > / < e^2 >
for i in range(1, n):
  if mms2[i-1] == 0 or mms2[i] == 0:
    dcr1 = 0
  else:
    dcr1 = .5 * (dbdus[i-1]*mms3[i-1]/mms2[i-1] + dbdus[i]*mms3[i]/mms2[i])
  crs1[i] = crs1[i-1] + .5*dcr1*4


# the second correction term
# - \int (log)' < e^2 >  = - log < e^2 >
for i in range(1, n):
  if mms2[i] > 0:
    if i < n/2:
      crs2[i] = -log(mms2[i]/mms2[n/2-1])
    elif i > n/2:
      crs2[i] = -log(mms2[i]/mms2[n/2+1])

#
# Corrections for < e exp(-bet/2 * e) > = 0
#

# recompute the density of states

sbeth[0] = log(2)
for i in range(1, n):
  ds = .5*4*(beth[i-1] + beth[i])
  if abs(ds) > 100: ds = 0
  sbeth[i] = sbeth[i-1] + ds

# the second correction term
# - \int (log)' < e^2 >*  = - log < e^2 >*
for i in range(1, n):
  if mmh2[i] > 0:
    if i < n/2:
      crh2[i] = -log(mmh2[i]/mmh2[n/2-1])
    elif i > n/2:
      crh2[i] = -log(mmh2[i]/mmh2[n/2+1])



# output 
for i in range(n):
  print "%g " * 23 % (
      ene[i], betx[i], mmx1[i], mmx2[i],
      lndos[i], sbetx[i], serrx[i], dbdux[i],
      crx1[i], crx2[i],  # 10,   $7 + $9 + $10 flat
      bets[i], mms2[i], mms3[i],
      sbets[i], serrs[i], dbdus[i],
      crs1[i], crs2[i],  # 18,  $15 + $17 + $18 flat
      beth[i], mmh2[i],
      sbeth[i], serrh[i],
      crh2[i],  # 23,  $22 + $23 flat!
      )



