#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 14:28:40 2020

Filenames for the ridgeline code.  Requires a ridges, problematic, distances,
CutOutCats, Ratios folders.  The .fits, .npys and data files can be obtained from
location; a cutouts folder for each of these is required.

@author: bonnybarkus
"""

path = 'Catalogue/DR1/'

# Load in files
#LofCat = str('Catalogue/DR2/data/13h_v0.3_new.fits')
LofCat = str('/Users/bonnybarkus/Documents/Ridgelines/New_Catalogue/Catalogue/data/all_eyeclass_4rms_v1.2_morph_agn_length_err.fits')
#CompCat = str('Catalogue/DR2/data/components-v0.3.fits')
CompCat = str('/Users/bonnybarkus/Documents/Ridgelines/New_Catalogue/Catalogue/data/LOFAR_HBA_T1_DR1_merge_ID_v1.2.comp.fits')
OptCat = str('Catalogue/data/pwerrunder1.fits')
OptCatdf = str('Catalogue/data/pwfull.txt')

# Data files
fitsfile = 'Catalogue/fits_cutouts/'
npyfile = 'Catalogue/rms4_cutouts/'
rms4 = 'Catalogue/rms4/'
rms4cutout = 'Catalogue/rms4_cutouts/'
fits = 'Catalogue/fits/'
fitscutout = 'Catalogue/fits_cutouts/'

# Ridgelines
TFC = path + 'TotalFluxCutFull.txt'
Probs = path + 'problematic/%s_image.png'
R1 = path + 'ridges/%s_ridge1.txt'
R2 = path + 'ridges/%s_ridge2.txt'
Rimage = path + 'ridges/%s_ridges%d.png'
psl = path + 'problematic_sources_list_Full.txt'

# SourceSearch
coc = path + 'CutOutCats/Cutout_Catalogue-%s.txt'
#Dists = 'Catalogue/DistancesFull/distances-%s.txt'
Position = path + 'Distances/Position_Info.txt'
RDists = path + 'Distances/Rdistances-%s.txt'
LDists = path + 'Distances/Ldistances-%s.txt'
NDist = path + 'Distances/NclosestLdistances-%s.txt'
NLLR = path + 'Ratios/NearestLofarLikelihoodRatios-%s.txt'
NRLR = path + 'Ratios/NearestRidgeLikelihoodRatios-%s.txt'
LLR = path + 'Ratios/LofarLikelihoodRatiosLR-%s.txt'
RLR = path + 'Ratios/RidgeLikelihoodRatiosLR-%s.txt'
NLRI = path + 'MagnitudeColour/Nearest30InfoLR-%s.txt'
LRI = path + 'MagnitudeColour/30InfoLR-%s.txt'
MagCO = path + 'MagCutOutCats/Cutout_Catalogue-%s.txt'

# Table Columns
LSN = 'Source_Name'  # Column of the LOFAR ID
LRA = 'RA'  # LOFAR catalogue position RA
LDEC = 'DEC'  # LOFAR catalogue position DEC
LredZ = 'z_best'
All = 'AllWISE'
OID = 'objID'
LMS = 'LM_size'
CSN = 'Source_Name'  # Component catalogue Source Name
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
PossHosts = 'Catalogue/PossibleHostsListDR1Full.csv'
#DR1Hosts = 'MagnitudeColourFull/HostMagnitude_Info.txt'
DR1HostsFull = path + 'MagnitudeColourFull/HostMagnitude_InfoFull.txt'
MCLR = path + 'MagnitudeColourFull/Nearest30AllLRW1band-%s.txt'
LR = path + 'MagnitudeColourFull/AllLRW1bandLR-%s.txt'
