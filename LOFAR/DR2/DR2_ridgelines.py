#!/usr/bin/env python 3
# coding: utf-8

# ## Ridgeline Code for DR2
# 
# This is a script version of notebook from BB containing the ridgeline code for DR2.  This is the preparation for running the likelihood ratio code and only needs to be run once for each dataset. 


# Imports

import time
import RidgelineFilesDR2 as RLF
import RLConstantsDR2 as RLC
from ridge_toolkitDR2 import AreaFluxes, ComponentsTable, CreateCutouts
from ridge_toolkitDR2 import DefineCutoutHDU, FindRidges, GetAvailableSources, GetCutoutArray
from ridge_toolkitDR2 import GetMaskedComp, TotalFluxSelector, TrialSeries
from os.path import exists
from astropy.table import Table
from warnings import simplefilter, resetwarnings
from sizeflux_tools import Flood

simplefilter('ignore') # there is a matplotlib issue with shading on the graphs

R = RLC.R
dphi = RLC.dphi

CompTable = ComponentsTable(RLF.CompCat)
ffo=Flood(CompTable)

# Run main ridgeline tasks to produce the directories of ridgelines

#if exists(RLF.TFC) == False:

print('Please wait output will show below.')
start_time = time.time()
available_sources=TotalFluxSelector(RLF.LofCat, ffo)
print('Time taken to calculate sizes = ' + str((time.time()-start_time)/(60*60)))
    
#available_sources = GetAvailableSources(str(RLF.TFC))
print('Number of sources in Sample = ',len(available_sources))
    
#CreateCutouts(available_sources)
#print('CutOuts created.')

print('Starting Ridgeline drawing process.')
start_time = time.time()
TrialSeries(available_sources, str(RLF.CompCat), R, dphi, ffo)
print('Time taken for Ridgelines to draw (Hours) = ' + str((time.time()-start_time)/(60*60)))
resetwarnings()
    
'''
else:
    print('Sample selected. Sizes and cutouts present. Ridgelines will now draw.')
    
    available_sources = GetAvailableSources(str(RLF.TFC))
    print('Number of sources in Sample = ' + str(available_sources.shape))
    
    print('Starting Ridgeline drawing process.')
    start_time = time.time()
    TrialSeries(available_sources, str(RLF.CompCat), R, dphi, CompTable)
    print('Time taken for Ridgelines to draw (Hours) = ' + str((time.time()-start_time)/(60*60)))
    resetwarnings()
'''
