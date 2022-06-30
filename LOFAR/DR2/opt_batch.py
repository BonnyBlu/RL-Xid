#!/usr/bin/python

'''
Filters radio and optical catalogues ready for processing, and batches radio catalogue into catalogues of 500 sources each in a separate sub-directory ready for DR2_setup.py.
'''

import numpy as np
from astropy.table import Table,vstack,unique
import sys
import os
import glob

root=sys.argv[1]
optfile="optical_filtered.fits"
ocat=Table.read(optfile)

dirs=glob.glob(root+"*/")

for d in dirs:
    rfile=d+d[:-1]+".fits"
    ofile=d+"optical.fits"
    rcat=Table.read(rfile)
    ramin=np.nanmin(rcat['RA'])
    ramax=np.nanmax(rcat['RA'])
    decmin=np.nanmin(rcat['DEC'])
    decmax=np.nanmax(rcat['DEC'])
    rl=ramin-0.5
    ru=ramax+0.5
    dl=decmin-0.5
    du=decmax+0.5
    subcat=ocat[(ocat['RA']>rl) & (ocat['DEC']>dl)]
    nsubcat=subcat[(subcat['RA']<ru) & (subcat['DEC']<du)]
    #filt &= dl<ocat['DEC']<du
    #newocat=ocat[filt]
    nsubcat.write(ofile)

