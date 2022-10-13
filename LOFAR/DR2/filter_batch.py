#!/usr/bin/python

'''
Filters radio and optical catalogues ready for processing, and batches radio catalogue into catalogues of 500 sources each in a separate sub-directory ready for DR2_setup.py.
'''

import numpy as np
from astropy.table import Table,vstack,unique
import sys
import os

rfil=sys.argv[1]
ofil=sys.argv[2]
outroot=sys.argv[3]
bpwd = '/data1/tap/ids/RL-Xid/P21_Test/'
try:
	os.mkdir(bpwd)
except:
	pass

rcat=Table.read(rfil)
ocat=Table.read(ofil)

# Filter radio catalogue for flux and size

rfcut=rcat[rcat['Total_flux']>10.0]
rmcut=rfcut[rfcut['Predicted_Width']>15.0]
rlcut=rfcut[rfcut['Predicted_Size']>15.0]
rsjoin=vstack([rmcut,rlcut])
rscut=unique(rsjoin)
print "Length of rm, rl: ",len(rmcut),len(rlcut)

lrad=len(rscut)
print "Length of filtered radio catalogue is "+str(lrad)

# Filter optical catalogue for WISE dets

ocut=ocat[ocat['UNWISE_OBJID']!="N/A"]

print "Length of filtered hosts catalogue is "+str(len(ocut))

# Batch radio catalogue

BNUM=25

nbatch=int(lrad/BNUM)+1

start=0
end=BNUM-1
for b in range(0,nbatch):
    bdir=outroot+str(b)
    bfolder = bpwd+bdir
    bname=bfolder+"/"+bdir+".fits"
    newrcat=rscut[start:end]
    #print bfolder, bdir, bname
    try:
        os.mkdir(bfolder)
    except:
        print "Directory "+bdir+" exists"
    newrcat.write(bname)
    start=end
    end=end+BNUM

# write full cats to parent directory

rscut.write("radio_filtered.fits")
ocut.write("optical_filtered.fits")


