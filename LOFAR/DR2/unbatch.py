 q#!/usr/bin/python

import numpy as np
from astropy.table import Table,vstack
import sys
import os
import glob

root=sys.argv[1]
outfile=root+"_allhosts.csv"
outfits=root+"_allhosts.fits"
outrand=root+"_randhosts.fits"

files=glob.glob(root+"*/hosts.csv")

lout=open(outfile,"w")
lout.write("Source_Name,RA,DEC,UNWISE_OBJID,UID_L,LRMagBoth,PossFail\n")
for d in files:
    lines=[l.rstrip() for l in open(d).readlines()]
    lines=lines[1:]
    #print lines +"\n"
    for l in lines:
        lout.write(l+"\n")
    #print d

lout.close()
    
# make random test sample of 300

intab=Table.read(outfile,format='ascii',delimiter=',')
gtab=intab[intab['LRMagBoth']>1.0]
nums=np.random.choice(len(gtab),300,replace=False)

keep=gtab[nums]

keep.write(outrand)

gtab.remove_columns(['UID_L','PossFail'])
gtab.rename_column('RA','optRA_RLC')
gtab.rename_column('DEC','optDEC_RLC')

gtab.write(outfits)

#print intab['LRMagBoth']
