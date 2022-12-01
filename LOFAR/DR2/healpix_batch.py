#!/usr/bin/python

from __future__ import print_function

'''
Filters radio and optical catalogues ready for processing, and batches radio catalogue into catalogues in a separate sub-directory ready for DR2_setup.py.

Input radio catalogue should be pre-filtered, e.g. all LGZ sources without optical IDs

This version by mjh -- filter using healpix rather than in order of catalogue entry

call with input radio catalogue (all suitable sources), optical catalogue (all IDs), component catalogue (all components)

'''

import numpy as np
from astropy.table import Table,vstack,unique
import sys
import os
from astropy_healpix import HEALPix
import astropy.units as u

rfil=sys.argv[1]
ofil=sys.argv[2]
cfil=sys.argv[3]
outroot='hp'
hp = HEALPix(nside=16)

print('Reading catalogues...')
rcat=Table.read(rfil)
ocat=Table.read(ofil)
ccat=Table.read(cfil)

# Filter radio catalogue for flux and size

rfcut=rcat[rcat['Total_flux']>10.0]
rmcutf=rfcut['LGZ_Width']>15.0
rlcutf=rfcut['LGZ_Size']>15.0
rscut=rfcut[rmcutf | rlcutf]
print("Length of rm, rl: ",np.sum(rmcutf),np.sum(rlcutf))

lrad=len(rscut)
print("Length of filtered radio catalogue is "+str(lrad))

hpix=hp.lonlat_to_healpix(rscut['RA']*u.deg,rscut['DEC']*u.deg)

uhp=sorted(np.unique(hpix))
print('There are',len(uhp),'healpixes in the radio catalogue')

# Filter optical catalogue for WISE dets

ocut=ocat[ocat['UNWISE_OBJID']!="N/A"]

print("Length of filtered hosts catalogue is "+str(len(ocut)))

ohpix=hp.lonlat_to_healpix(ocat['RA']*u.deg,ocat['DEC']*u.deg)

chpix=hp.lonlat_to_healpix(ccat['RA']*u.deg,ccat['DEC']*u.deg)

for i,pix in enumerate(uhp):
    print('Doing healpix',i,pix)
    bdir=os.path.join(os.getenv('ID_RESULTS'),outroot+'_'+str(i))
    bname=bdir+'/'+outroot+'_'+str(pix)+'.fits'
    #print("Bdir",bdir)
    #print("outrout",outroot)
    #print("pix",pix)
    #print("Bname",bname)
    newrcat=rscut[hpix==pix]
    try:
        os.mkdir(bdir)
    except:
        print("Directory "+bdir+" exists")
    newrcat.write(bname,overwrite=True)
    # the optical and components catalogues are for this pix and its
    # immediate neighbours -- wasteful but ensures no boundary effects
    neighbours=np.ndarray.flatten(hp.neighbours(pix))
    cmask=np.zeros_like(chpix,dtype=bool)
    omask=np.zeros_like(ohpix,dtype=bool)
    for npix in list(neighbours)+[pix]:
        cmask|=(chpix==npix)
        omask|=(ohpix==npix)
    newccat=ccat[cmask]
    newccat.write(bdir+'/'+'components.fits',overwrite=True)
    ocat[omask].write(bdir+'/'+'optical.fits',overwrite=True)

rscut.write(os.path.join(os.getenv('ID_RESULTS'),"radio_filtered.fits"),overwrite=True)
#ocut.write("optical_filtered.fits",overwrite=True)


