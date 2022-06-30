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

sourcecat = sys.argv[1]
compcat = sys.argv[2]

sources=Table.read(sourcecat)
comps=Table.read(compcat)

for row in sources:
    ssource=row['Source_Name']
    ssource=ssource.rstrip()
    row['Source_Name']=ssource

sources.write('radio.fits',overwrite=True)

for nrow in comps:
    cname=nrow['Component_Name']
    pname=nrow['Parent_Source']
    cname=cname.rstrip()
    pname=pname.rstrip()
    nrow['Component_Name']=cname
    nrow['Parent_Source']=pname

comps.write('../components.fits',overwrite=True)
