#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 12:46:17 2019

Constants for the ridgeline toolkit, Host_Distances and SourceSearch.

@author: bonnybarkus
"""

##  Allows printstatement to show progress if set to 'True'
debug = 'True'

##  AreaFluxes
KerW = 4  ##  Weight of the kernel in the convolution

##  ConfigurePlots
lw = 1  ##  Linewidth of graphs
col = 'green'  ##  line.color has to be a string
ms = 3  ##  markersize
mt = 'serif'  ##  mathtext  has to be a string
ff = 'serif'  ##  font.family has to be a string
fs = 12  ##  fontsize
ts = 15  ##  titlesize
ls = 15  ##  labelsize
ym = 0.5  ##  ymargin
figs = 7  ##  figsize draws a square if insert twice
fonts = 'large'  ##  legend fontsize has to be a string
sav = 'eps'  ##  save format has to be a string
pcol = 'auto'  ## shading type

##  ErodedMaxima
TRel = 0.5  ##  Value Relative threshold must be greater than in the erosion.
pThresh = 0.2  ##  Peak threshold value for Erosion
mindis = 5  ##  Minimum dist between peaks in Erosion
Octm = 3  ##  The m parameter of the Octagon function used to erode
Octn = 3  ##  The n parameter of the octagon function used to erode

##  FindRidges
Lcone = 85  ##  The angle in degrees of the larger cone
MaxLen = 0.95 ##  Multiplier of Source size to determine max RL length
ipit = 6  ##  Max number of iterations in the initial point search for directions

##  FloodFill, GetAvailableSources, CreateCutOutCat
rdel = -0.000305556  ## The equivalent pixel value for RA in FITS file in degrees (DON'T CHANGE)
ddel = 0.000305556  ## The equivalent pixel value for DEC in FITS file in degrees (DON'T CHANGE)

##  GetRidgePoint
Jlim = 40/95  ##  Mulitpler of RL Length Max that it can jump under
MaxJ = 1  ##  Multiplier of source size the RL can jump to
JMax = 1  ##  Multiplier of max RL len (for RL total to be subtracted from)

##  MaskedCentre
Rad = 2.5 ##  Radius of the circle to mask the centre

##  TotalFluxSelector
LGZA = 1.  ##  Number of Components to filter by if active.  Needs dots.

##  TrialSeries, CreateCutouts
Start = 0  ##  The index to start the selection of sources from
Step = 1  ##  The fraction step size of the function if 1 runs all, 2 50%, 4 25% etc
ImSize = 1  ##  The fraction of the source the final image is cut down to

##  GetRMSArray
nSig = 4.0  ##  The multiple of sigma for the RMS reduction

##  HostRFromLofar
SigAst = 0.6 ##  The astrometric error, changing between a radio and optical catalogue (arcsec)

##  Lambda
optcount = 26569141 ##  The number of optical sources in the catalogue
catarea =  11431476680.687202  ##  The area of the catalogue in arcseconds squared (424 x 3600 x 3600 for LOFAR)

##  Cutout Creation
rsize = 1/60  #  The distance away along the RA in degrees to form the sub-catalogue
dsize = 1/120  #  The distance away along the DEC in degrees to form the sub-catalogue

##  LikelihoodRatios
radraerr = 0.2  ## The RA radio on the ridgeline (defined as zero until I know better)
raddecerr = 0.2  ## The DEC radio error on the ridgeline
UniWErr = 0.2
UniLErr = 0.1
LRAd = 160.0 # float, degrees, RA down the lower value of RA of the Optical sky area
LRAu = 232.0 # float, degrees, RA up the upper value of RA of the Optical sky area
LDECd = 42.0 # float, degrees, DEC down the lower value of DEC of the Optical sky area
LDECu = 62.0 # float, degrees, DEC up the upper value of the DEC of the Optical sky area
bw = 0.2  # Bandwidth of KDE
SigLC = 0.2 # Lofar centre distribution
Lth = 0.0 # Threshold for declaring a PossFail