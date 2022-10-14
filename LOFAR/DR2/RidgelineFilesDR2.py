#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 14:28:40 2020

Filenames for the ridgeline code. This is a template that is completed by the DR2_setup.py script

@author: bonnybarkus ; mangled a bit by JHC
"""
import os

# Load in files -- added by setup script

indir=os.getenv('RLDIR')

LofCat = indir+'/radio.fits'
CompCat = indir+'/components.fits' # changed by MJH -- each dir has
                                   # its own components file now
OptCat = indir+'/optical.fits'
OptCatdf = indir+'/optical.txt'
PossHosts = indir+'/hosts.csv'
tmpdir = indir+'/'

# Data files
fitsfile = indir+'/fits_cutouts/'
npyfile = indir+'/rms4_cutouts/'
rms4 = indir+'/rms4/'
rms4cutout = indir+'/rms4_cutouts/'
fits = indir+'/fits/'
fitscutout = indir+'/fits_cutouts/'

# Ridgelines
TFC = indir+'/total_flux_cutWorkingSet.txt'
Probs = indir+'/problematic/%s_image.png'
R1 = indir+'/ridges/%s_ridge1.txt'
R2 = indir+'/ridges/%s_ridge2.txt'
Rimage = indir+'/ridges/%s_ridges%d.png'
psl = indir+'/problematic/problematic_sources_listWorkingSet.txt'

# SourceSearch
coc = indir+'/CutOutCats/Cutout_Catalogue-%s.txt'
#Dists = 'Catalogue/DistancesFull/distances-%s.txt'
Position = indir+'/Distances/Position_Info.txt'
RDists = indir+'/Distances/Rdistances-%s.txt'
LDists = indir+'/Distances/Ldistances-%s.txt'
NDist = indir+'/Distances/NclosestLdistances-%s.txt'
NLLR = indir+'/Ratios/NearestLofarLikelihoodRatios-%s.txt'
NRLR = indir+'/Ratios/NearestRidgeLikelihoodRatios-%s.txt'
LLR = indir+'/Ratios/LofarLikelihoodRatiosLR-%s.txt'
RLR = indir+'/Ratios/RidgeLikelihoodRatiosLR-%s.txt'
NLRI = indir+'/MagnitudeColour/Nearest30InfoLR-%s.txt'
LRI = indir+'/MagnitudeColour/30InfoLR-%s.txt'
MagCO = indir+'/MagCutOutCats/Cutout_Catalogue-%s.txt'

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
OptMagA = 'MAG_W1'  # Magnitude from AllWISE
MagAErr = ''
OptMagP = 'MAG_R'  # Magnitude from PanSTARRS
MagPErr = ''
PossRA = 'RA'  # Possible Optical counterpart RA
PossDEC = 'DEC'  # Possible Optical counterpart DEC
PRAErr = 0.2  # Error on possible Optical counterpart RA
PDECErr = 0.2  # Error on possible Optical counterpart DEC
IDW = 'UNWISE_OBJID' # WISE ID
IDP = 'UID_L' # PS ID - Unique ID Legacy
ID3 = 'UNWISE_OBJID' # ID to be taken one or the other of WISE or Legacy
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
#Odata = '/beegfs/lofar/jcroston/surveys/dr2_hosts/pwfull.txt' # Original DR1 optical txt file [not needed?]
DR1Hosts = '/beegfs/lofar/jcroston/surveys/dr2_hosts/testing/HostMagnitude_Info.txt'
DR1HostsFull = '/beegfs/lofar/jcroston/surveys/dr2_hosts/testing/HostMagnitude_InfoFull.txt'
MCLR = indir+'/MagnitudeColour/Nearest30AllLRW1band-%s.txt'
LR = indir+'/MagnitudeColour/AllLRW1bandLR-%s.txt'

