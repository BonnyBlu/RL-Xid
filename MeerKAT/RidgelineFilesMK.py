#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 14:28:40 2020

Filenames for the ridgeline code.  Requires a ridges, problematic, distances,
CutOutCats, Ratios folders.  The .fits, .npys and data files can be obtained from
location; a cutouts folder for each of these is required.

@author: bonnybarkus
"""

#path = '/home/bbarkus/MeerKAT/'

# Load in files
#LofCat = str('Catalogue/DR2/data/13h_v0.3_new.fits')
LofCat = str('/home/bbarkus/MeerKAT/xmmlss_sources_v1.fits')
#CompCat = str('Catalogue/DR2/data/components-v0.3.fits')
CompCat = str('/home/bbarkus/MeerKAT/xmmlss_components_v1.fits')
OptCat = str('Catalogue/data/pwerrunder1.fits')
OptCatdf = str('Catalogue/data/pwfull.txt')

# Data files
fitsfile = 'fits_cutouts/'
npyfile = 'rms4_cutouts/'
rms4 = 'rms4/'
rms4cutout = 'rms4_cutouts/'
fits = 'fits/'
fitscutout = 'fits_cutouts/'

# Ridgelines
TFC = 'TotalFluxCutMK.txt'
Probs = 'problematicMK/%s_image.png'
R1 = 'ridgesMK/%s_ridge1.txt'
R2 = 'ridgesMK/%s_ridge2.txt'
Rimage = 'ridgesMK/%s_ridges%d.png'
psl = 'problematic_sources_listMK.txt'

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
SSN = 'Source_Name'  # Source catalogue Source Name
STF = 'Total_flux'  # Source catalogue total flux column
SRA = 'RA'  # Source catalogue position RA
SRAE = 'E_RA'  # Source catalogue RA error
SDEC = 'DEC'  # Source catalogue position DEC
SDECE = 'E_DEC'  # Source catalogue DEC error
SASS = 'Assoc'  # Source catalogue number of component associations column
CSN = 'Parent_Source'  # Component catalogue Source Name
CCN = 'Component_Name'  # Component catalogue Component Name
CTF = 'Total_flux'  # Component catalogue total flux column
CRA = 'RA'  # Component catalogue position RA
CDEC = 'DEC'  # Component catalogue position DEC
CMAJ = 'Maj'  # Component Major axis column
CMIN = 'Min'  # Component Minor axis column
CPA = 'PA'  # Component rotational angle column
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
OptMagA = 'W1mag'  # Magnitude from AllWISE
MagAErr = 'W1magErr'
OptMagP = 'i'  # Magnitude from PanSTARRS
MagPErr = 'iErr'

# Magnitude and Colour Likelihood Ratio
Odata = 'data/pwfull.txt' # Original DR1 optical txt file
#MLR = 'Catalogue/MagnitudeColourFull/LikelihoodRatios-%s.txt'
PossHosts = 'Catalogue/PossibleHostsList.csv'
DR1Hosts = 'MagnitudeColour/HostMagnitude_Info.txt'
DR1HostsFull = 'MagnitudeColour/HostMagnitude_InfoFull.txt'
MCLR = 'MagnitudeColour/Nearest30AllLRW1band-%s.txt'
LR = 'MagnitudeColour/AllLRW1bandLR-%s.txt'
