#!/usr/bin/python

from __future__ import print_function

'''
Filters radio and optical catalogues ready for processing, and batches radio catalogue into catalogues of 500 sources each in a separate sub-directory ready for DR2_setup.py.

Input radio catalogue should be pre-filtered, e.g. all LGZ sources without optical IDs
'''

import numpy as np
from astropy.table import Table,vstack,unique
import sys
import os

rfil=sys.argv[1]
ofil=sys.argv[2]
outroot=sys.argv[3]

rcat=Table.read(rfil)
ocat=Table.read(ofil)

# Filter radio catalogue for flux and size

rfcut=rcat[rcat['Total_flux']>10.0]
rmcutf=rfcut['Maj']>15.0
rlcutf=rfcut['LGZ_Size']>15.0
rscut=rfcut[rmcutf | rlcutf]
print("Length of rm, rl: ",np.sum(rmcutf),np.sum(rlcutf))

lrad=len(rscut)
print("Length of filtered radio catalogue is "+str(lrad))

# Filter optical catalogue for WISE dets

ocut=ocat[ocat['UNWISE_OBJID']!="N/A"]

print("Length of filtered hosts catalogue is "+str(len(ocut)))

# Batch radio catalogue

BNUM=500

nbatch=int(lrad/500.0)+1

start=0
end=499
for b in range(0,nbatch):
    bdir=outroot+str(b)
    bname=bdir+"/"+bdir+".fits"
    newrcat=rscut[start:end]
    try:
        os.mkdir(bdir)
    except:
        print("Directory "+bdir+" exists")
    newrcat.write(bname,overwrite=True)
    start=end
    end=end+500

# write full cats to parent directory

rscut.write("radio_filtered.fits",overwrite=True)
ocut.write("optical_filtered.fits",overwrite=True)


