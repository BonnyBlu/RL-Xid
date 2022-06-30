#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 14:28:40 2020

Filenames for the ridgeline code. This is a template that is completed by the DR2_setup.py script

@author: bonnybarkus ; mangled a bit by JHC
"""

# Load in files -- added by setup script

# Data files
fitsfile = 'fits_cutouts/'
npyfile = 'rms4_cutouts/'
rms4 = 'rms4/'
rms4cutout = 'rms4_cutouts/'
fits = 'fits/'
fitscutout = 'fits_cutouts/'

# Ridgelines
TFC = 'total_flux_cutWorkingSet.txt'
Probs = 'problematic/%s_image.png'
R1 = 'ridges/%s_ridge1.txt'
R2 = 'ridges/%s_ridge2.txt'
Rimage = 'ridges/%s_ridges%d.png'
psl = 'problematic/problematic_sources_listWorkingSet.txt'

# SourceSearch
coc = 'CutOutCats/Cutout_Catalogue-%s.txt'
#Dists = 'Catalogue/DistancesFull/distances-%s.txt'
Position = 'Distances/Position_Info.txt'
RDists = 'Distances/Rdistances-%s.txt'
LDists = 'Distances/Ldistances-%s.txt'
NDist = 'Distances/NclosestLdistances-%s.txt'
NLLR = 'Ratios/NearestLofarLikelihoodRatios-%s.txt'
NRLR = 'Ratios/NearestRidgeLikelihoodRatios-%s.txt'
LLR = 'Ratios/LofarLikelihoodRatiosLR-%s.txt'
RLR = 'Ratios/RidgeLikelihoodRatiosLR-%s.txt'
NLRI = 'MagnitudeColour/Nearest30InfoLR-%s.txt'
LRI = 'MagnitudeColour/30InfoLR-%s.txt'
MagCO = 'MagCutOutCats/Cutout_Catalogue-%s.txt'

# Table Columns
OptMagA = 'W1mag'  # Magnitude from AllWISE
MagAErr = 'W1magErr'
OptMagP = 'i'  # Magnitude from PanSTARRS
MagPErr = 'iErr'
PossRA = 'ra'  # Possible Optical counterpart RA
PossDEC = 'dec'  # Possible Optical counterpart DEC
PRAErr = 'raErr'  # Error on possible Optical counterpart RA
PDECErr = 'decErr'  # Error on possible Optical counterpart DEC
IDW = 'AllWISE' # WISE ID
IDP = 'objID' # PS ID - Unique ID Legacy
#ID3 = 'ID' # ID to be taken one or the other of WISE or Legacy
LRMA = 'LRMagA'
LRMP = 'LRMagR'
LRMC = 'LRMagBoth'
LSN = 'Source_Name'  # Column of the LOFAR ID
LRA = 'RA'  # LOFAR catalogue position RA
LDEC = 'DEC'  # LOFAR catalogue position DEC
LredZ = 'z_best'
All = 'AllWISE'
OID = 'objID'
LMS = 'LM_size'

# Magnitude and Colour Likelihood Ratio
Odata = '/beegfs/lofar/jcroston/surveys/dr2_hosts/pwfull.txt' # Original DR1 optical txt file [not needed?]
DR1Hosts = '/beegfs/lofar/jcroston/surveys/dr2_hosts/HostMagnitude_Info.txt'
DR1HostsFull = '/beegfs/lofar/jcroston/surveys/dr2_hosts/HostMagnitude_InfoFull.txt'
MCLR = 'MagnitudeColour/Nearest30AllLRW1band-%s.txt'
LR = 'MagnitudeColour/AllLRW1bandLR-%s.txt'
LofCat = spring60-65_bb.fitsCompCat = components-v0.1.fitsOptCat = optical_wisedets2.fitsPossHosts = spring60-65_bb_RLhosts.csv