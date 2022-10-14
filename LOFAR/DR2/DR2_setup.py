#!/usr/bin/env python 3
# coding: utf-8

# ## Ridgeline Code for DR2
# 
# This script sets up the input files and directory structure to run the ridgeline and LR steps on a particular fits catalogue with corresponding components cat and optical potential hosts catalogue.
#
# It requires IMAGEDIR to be set and the lotss-catalogue utils to be on its path. It should be run in a new directory containing the input catalogues, and will make a set of directories for the RLC intermediate products. It will also write out an updated RidgelineFilesDR2.py with the correct input and output names for the host ID step.

# Imports

from __future__ import print_function

import numpy as np
import sys
import os
from astropy.table import Table
from astropy.io import fits
import glob
from subim import extract_subim
from overlay import find_noise_area,find_noise
from download_image_files import LofarMaps
from astropy.coordinates import SkyCoord
import astropy.units as u
import time
import RidgelineFilesDR2 as RLF
import RLConstantsDR2 as RLC
from os.path import exists
from warnings import simplefilter, resetwarnings
simplefilter('ignore') # there is a matplotlib issue with shading on the graphs

def get_fits(fra,fdec,fsource,fsize):

    sc=SkyCoord(fra*u.deg,fdec*u.deg,frame='icrs')
    s=sc.to_string(style='hmsdms',sep='',precision=2)
    name=fsource
    newsize=2.5*fsize/3600.0

    lm=LofarMaps()
    mosname=lm.find(fra,fdec)
    filename=os.environ['IMAGEDIR']+'/'+mosname
    hdu=extract_subim(filename,fra,fdec,newsize)
    if hdu is not None:
        hdu.writeto('fits/'+name+'.fits',overwrite=True)
        flag=0
    else:
        print('Cutout failed for',fsource)
        flag=1

    return flag

sourcecat = sys.argv[1]
compcat = sys.argv[2]
#compcat = sys.argv[2]
#hostcat = sys.argv[3]
#outroot = sys.argv[4]

#inridge="RidgelineFiles_temp.py"
newdirs=['fits','rms4','fits_cutouts','rms4_cutouts','Distances','MagnitudeColour','Ratios','CutOutCats','MagCutOutCats','badsources_output','ridges','problematic','cutouts']

path=os.getcwd()

for d in newdirs:
    newd=path+'/'+d
    try:
        os.mkdir(newd)
    except:
        # dir already exists
        print("Directory",newd,"already exists, cleaning it out")
        os.system("rm "+newd+"/*")
    else:
        print("Made directory ",newd)

sources=Table.read(sourcecat)
comps=Table.read(compcat)
# Extract cutouts and thresholded npy arrays

for row in sources:
    maj=row['Predicted_Width']
    lgzsiz=row['LGZ_Size']
    ssource=row['Source_Name']
    ssource=ssource.rstrip()
    row['Source_Name']=ssource
    sra=row['RA']
    sdec=row['DEC']
    #optra=row['optRA']
    #optdec=row['optDec']
    flux=row['Peak_flux']
    rms=row['Isl_rms']
    
    if np.isnan(maj):
        ssize=lgzsiz
    else:
        ssize=maj
   
    if ssize<20.0:
        ssize=20.0

    flag=get_fits(sra,sdec,ssource,ssize)

    cutout=path+'/fits/'+ssource+'.fits'
    lmd=LofarMaps()
    dname=os.environ['IMAGEDIR']+'/'+lmd.find(sra,sdec)
    if os.path.isfile(dname):
        lhdu=fits.open(dname)
        if flag==0:
            nlhdu=fits.open(cutout) 
            d=nlhdu[0].data
            peak=flux
            dyncut=50.0
            ratio=peak/dyncut
            print(peak,rms,ratio)
            if rms<ratio:
                print("rms < ratio")
                thres=(1e-3)*ratio
                thres2=ratio
                thres3=ratio
            else:
                print("rms > ratio")
                print(rms)
                thres=(1e-3)*4.0*rms

            d[d<thres]=np.nan
            mtest=np.nanmax(d)
            print("Max val of thresholded array is:",mtest)
            np.save(path+"/rms4/"+ssource+'.npy',d)
 
print("Completed generating fits and thresholded npy cutouts.")

sources.write('radio.fits',overwrite=True)

for nrow in comps:
    cname=nrow['Component_Name']
    pname=nrow['Source_Name']
    cname=cname.rstrip()
    pname=pname.rstrip()
    nrow['Component_Name']=cname
    nrow['Source_Name']=pname

#comps.write('../components.fits')

# Append input and output lines to RidgelineFiles template
'''
rlines=[l.rstrip().split(",") for l in open(inridge).readlines()]

rfile=open(inridge,"a")

rfile.write("LofCat = \""+sourcecat+"\"\n")
rfile.write("CompCat = \""+compcat+"\"\n")
rfile.write("OptCat = \""+hostcat+"\"\n")
rfile.write("PossHosts = \""+outroot+"_RLhosts.csv\"\n")

rfile.close()

cpcmd="cp "+sourcecat+" radio.fits"
os.system(cpcmd)
'''
