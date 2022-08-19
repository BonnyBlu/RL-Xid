#!/usr/bin/env python3

# coding: utf-8

# ## Likelihood Ratio for DR2 ##
# 
# This notebook has the running of the Likelihood Ratio Code on DR2.  Once the output from the Ridgeline Code has been produced this code can be used to determine the list of possible hosts for the DR2 data set.

# In[ ]:


# Imports
import numpy as np
from astropy.io import fits
import os
from astropy.table import Table
import pandas as pd
from os.path import exists
from sklearn.neighbors import KernelDensity
import RidgelineFilesDR2 as RLF
import RLConstantsDR2 as RLC
from ridge_toolkitDR2 import DefineCutoutHDU, GetAvailableSources, GetCutoutArray
from SourceSearchDR2 import GetSourceList, SourceInfo, filter_sourcelist, CreateLDistTable, NClosestDistances, GetCatArea, TableFromLofar, CreateSubCat, CreateCutOutCat, LikelihoodRatios
import sys

fitsfile=RLF.TFC.replace('.txt','.fits')
available_sources=Table.read(fitsfile)

# Load in the optical/IR and LOFAR catalogues, form tables and df and save as text files 
# Check columns esp on the Opt/IR
#OptTable = TableOfSources(str(RLF.OptCat))
OptTable = Table.read(RLF.OptCat)
# There aren't zillions of columns so fine to use whole thing
LofarTable = TableFromLofar(RLF.LofCat)
Lofardf = LofarTable.to_pandas()
Optdf = OptTable.to_pandas()
#Optdf.to_csv(RLF.OptCatdf, columns = [str(RLF.IDW), str(RLF.IDP), str(RLF.PossRA), str(RLF.PossDEC), str(RLF.OptMagA), str(RLF.OptMagP)], header = True, index = False)



# Set up source_list to be used in all of the following fuctions/cells
probfile = RLF.psl
source_list = GetSourceList(available_sources, probfile)

# source_list is now a filtered table

zerococ=[]
for source in source_list:
    source_name = source['Source_Name']
    lofarra = source['RA']
    lofardec = source['DEC']
    sizepix = source['Size']

    size = sizepix * RLC.ddel # convert size in pixels to degrees
    subcat = Optdf[(np.abs(Optdf[str(RLF.PossRA)] - lofarra) * np.cos(lofardec * np.pi / 180.0) < size) & (np.abs(Optdf[str(RLF.PossDEC)] - lofardec) < size)].copy()

    # Insert the uniform optical position error if required             

    subcat['raErr'] = np.where(np.isnan(subcat[str(RLF.OptMagP)]), RLC.UniWErr, RLC.UniLErr)
    subcat['decErr'] = np.where(np.isnan(subcat[str(RLF.OptMagP)]), RLC.UniWErr, RLC.UniLErr)

    subcat.to_csv(RLF.MagCO %source_name, columns = [str(RLF.IDW), str(RLF.IDP), str(RLF.PossRA), str(RLF.PossDEC), 'raErr', 'decErr', str(RLF.OptMagA), str(RLF.OptMagP)], header = True, index = False)
    print("Subcat written for source",source_name)

            

# Looping through all successful sources to create the cutoutcat .txt files.  The distance away to form
# the sub-catalogue is set in RLConstants and is currently set to 1 arcmin RA and 0.5 arcmin DEC.
# Only needs to be done once

source_count = 0 ## Keeps track of where the loop is
#print("Length of source list before COC loop: ",len(source_list))

for source in source_list:
    source_name=source['Source_Name']
    size = source['Size']
    #lofar_ra, lofar_dec = SourceInfo(source, LofarTable)[:2] # why was this needed?
    lofar_ra = source['RA']
    lofar_dec = source['DEC']

    subcat1 = CreateSubCat(OptTable, lofar_ra, lofar_dec)
    print('Optical subcat length is',len(subcat1))
    
    # Insert the uniform optical position error if required
    subcatdf = subcat1.to_pandas()

    # Insert the uniform optical position error if required             

    subcatdf['raErr'] = np.where(np.isnan(subcatdf[str(RLF.OptMagP)]), RLC.UniWErr, RLC.UniLErr)
    subcatdf['decErr'] = np.where(np.isnan(subcatdf[str(RLF.OptMagP)]), RLC.UniWErr, RLC.UniLErr)
    subcat2 = Table.from_pandas(subcatdf)
    #print("About to get COC for ",source)
    cutoutcat,lcat = CreateCutOutCat(source_name, LofarTable, subcat2, lofar_ra, lofar_dec, size)
    source_count += 1
    print("Source Number = ",source_count)
    if lcat==0:
        print("Removing source for zero-length cutout cat: ",source_name)
        zerococ.append(source_name)

source_list=filter_sourcelist(source_list,zerococ)
        
print('Number of viable sources = ',len(source_list))

if len(source_list)==0:
    sys.exit(0)

# Create a table of R distance information from LOFAR Catalogue position
# Only needs to be done once - need to exclude probs?
for source in source_list:
    print("Creating distance table for:",source['Source_Name'])
    CreateLDistTable(source['Source_Name'],source_list)


# Find the 30 closest sources for each ridgeline
# Only needs to be done once
n = 30
print('Finding closest distances')
NClosestDistances(source_list, LofarTable, n)

# In[ ]:

# get min and max RA and dec directly from radio catalogue

area=GetCatArea(RLF.OptCat)
print('Catalogue area is',area)

# Generating the likelihood ratios for all possible close sources, for each drawn
# ridgeline, using the R distance from LOFAR Catalogue position and the ridgeline.
# Only needs to be done once
print('Finding LRs')
LikelihoodRatios(source_list)

# CREATE NEAREST 30 INFO
# Load in the three text files for each source, join all the table information together and save
# Only needs to be done once

for asource in source_list:
    source=asource['Source_Name']
    print('Creating nearest 30 tables for',source)

    LofarLR = pd.read_csv(RLF.LLR %source, header = 0)
    RidgeLR = pd.read_csv(RLF.RLR %source, header = 0, usecols = ['Ridge_LR'])
    MagCutOut = pd.read_csv(RLF.MagCO %source, header = 0, usecols = [RLF.IDW, RLF.IDP, RLF.PossRA, RLF.OptMagA, RLF.OptMagP])
    MagCutOut[str(RLF.PossRA)] = MagCutOut[str(RLF.PossRA)].apply(lambda x: round(x, 7))
            
    All_LR = LofarLR.join(RidgeLR['Ridge_LR'])
    # Changed combined to use just the lofar value if the ridge value is nan
    All_LR['Multi_LR'] = np.where(~np.isnan(All_LR['Ridge_LR']), All_LR['Lofar_LR'].astype(np.float64).multiply(All_LR['Ridge_LR'].astype(np.float64), axis = 'index'), All_LR['Lofar_LR'].astype(np.float64))
        
    All_LR.columns=['LofarRDis', 'Lofar_LR', str(RLF.PossRA), str(RLF.PossDEC), 'Ridge_LR', 'Multi_LR']
    All_LR[str(RLF.PossRA)] = All_LR[str(RLF.PossRA)].apply(lambda x: round(x, 7))
            
    MagLR = All_LR.merge(MagCutOut, on = str(RLF.PossRA))
            
    MagLR.to_csv(RLF.LRI %source, columns = ['LofarRDis', 'Lofar_LR', str(RLF.PossRA), str(RLF.PossDEC), 'Ridge_LR', 'Multi_LR', str(RLF.IDW), str(RLF.IDP), str(RLF.OptMagP), str(RLF.OptMagA)], header = True, index = False)


print('Colour LR calculations')


# Load in parent population host info
Hosts = pd.read_csv(str(RLF.DR1Hosts), usecols = [ 'Source_Name', 'AllWISE', 'Host_RA', 'Host_DEC', 'W1mag', 'i'], header = 0)
Hosts['r'] = Hosts['i'].apply(lambda y: y + 0.33)
Hosts['Colour'] = Hosts['r'].astype(np.float64).subtract(Hosts['W1mag'].astype(np.float64), axis = 'index')


# In[47]:


# Create the colour column and sample the df
#pwfulldf = pd.read_csv(RLF.OptCatdf, header = 0, usecols = [RLF.IDW,RLF.OptMagP,RLF.OptMagA])
pwfulldf=Optdf # already in memory
# Temporary update to keep in WISE-only galaxies
ColourPW = pwfulldf[~np.isnan(pwfulldf[RLF.OptMagP]) & ~np.isnan(pwfulldf[RLF.OptMagA])].copy()
ColourPW.reset_index(drop = True, inplace = True)
ColourPW['Colour'] = ColourPW[RLF.OptMagP].subtract(ColourPW[RLF.OptMagA], axis = 'index')
ColSam = ColourPW.sample(50000, replace = False)



'''
# Skyarea Covered byt he LOFAR data set
area = (np.deg2rad(RLC.LRAu) - np.deg2rad(RLC.LRAd)) * (np.sin(np.deg2rad(RLC.LDECu)) - np.sin(np.deg2rad(RLC.LDECd))) * np.rad2deg(3600)**2
'''


# Not running on i-band so only need the W1 band cells
hh, ww1 = np.mgrid[Hosts['Colour'].min() : Hosts['Colour'].max() : 0.05, Hosts['W1mag'].min() : Hosts['W1mag'].max() : 0.05]
h_sample = np.vstack([ww1.ravel(), hh.ravel()]).T
h_train = np.vstack([Hosts['W1mag'], Hosts['Colour']]).T
kde_h = KernelDensity(kernel = 'gaussian', bandwidth = RLC.bw)
kde_h.fit(h_train)
prob_h = np.exp(kde_h.score_samples(h_sample))
norm_h = len(Hosts['W1mag'])/np.sum(prob_h)


# In[50]:


oo, ww2 = np.mgrid[ColSam['Colour'].min() : ColSam['Colour'].max() : 0.05, ColSam[RLF.OptMagA].min() : ColSam[RLF.OptMagA].max() : 0.05]
o_sample = np.vstack([ww2.ravel(), oo.ravel()]).T
o_train = np.vstack([ColSam[RLF.OptMagA], ColSam['Colour']]).T
kde_o = KernelDensity(kernel = 'gaussian', bandwidth = RLC.bw)
kde_o.fit(o_train)
prob_o = np.exp(kde_o.score_samples(o_sample))
norm_o = len(ColSam[RLF.OptMagA])/np.sum(prob_o)


def GetLR(fr, qm, nm):
    lr = (fr * qm) / nm
    return lr

def GetLR2(fr, qm, nm,debug=False):
    if debug: print("fr, qm and nm are: ",fr,qm,nm)
    # we want density in units of per typical source length
    # squared for LOFAR f(R) and sqrt(density) in arcsec
    #for RL f(R) - take 60 arcsec
    lr = (fr * qm) / ((nm**1.5)/(np.sqrt(area)*area/(60.0*60.0)))
    return lr

def Getqmc(m, c):
    qmc = np.exp(kde_h.score_samples(np.array([m, c]).reshape(1, -1)))
    return qmc * norm_h

def Getnmc(m, c):
    nmc = np.exp(kde_o.score_samples(np.array([m, c]).reshape(1, -1)))
    return nmc * norm_o


# Calculating the LR from the text files for the W1 band hosts

# Removed the line where you deal with taking the W1 value if the r band value was 0.

# try just fixing the colours here

newdrop=[]
for asource in source_list:
    source=asource['Source_Name']
            
    #MLRhw = pd.read_csv(str(RLF.NLRI) %source, header = 0, usecols = ['AllWise', 'LofarRDis', 'ra', 'dec', 'Lofar_LR', 'Ridge_LR', 'SBLR', 'Multi_LR',  'i', 'W1mag'])
    MLR = pd.read_csv(str(RLF.LRI) %source, header = 0, usecols = ['LofarRDis', str(RLF.PossRA), str(RLF.PossDEC), str(RLF.IDW), str(RLF.IDP), str(RLF.OptMagP), str(RLF.OptMagA), 'Multi_LR'])
    MLR['Col'] = MLR[RLF.OptMagP].subtract(MLR[RLF.OptMagA], axis = 'index')
    MLR['Colour'] = MLR.apply(lambda row: np.where(row[RLF.OptMagP]>0.0,row[RLF.OptMagP]-row[RLF.OptMagA],RLC.meancol),axis=1).astype(np.float128)
    MCLT = MLR[~np.isnan(MLR['Colour'])].copy()
    MCLR = MCLT[~np.isnan(MLR[RLF.OptMagA])].copy()
    print("Source is",source,"and MCLR has length",len(MCLR))
    #MCLR['MCLLR'] = MCLR.apply(lambda row: GetLR(row['Lofar_LR'], Getqmc(row['W1mag'], row['Colour']), Getnmc(row['W1mag'], row['Colour'])), axis = 1).astype(np.float128)
    #MCLR['MCRLR'] = MCLR.apply(lambda row: GetLR(row['Ridge_LR'], Getqmc(row['W1mag'], row['Colour']), Getnmc(row['W1mag'], row['Colour'])), axis = 1).astype(np.float128)
    #MCLRhw['MCSBLR'] = MCLRhw.apply(lambda row: GetLR(row['SBLR'], Getqmcw(row['W1mag'], row['Colour']), Getnmcw(row['W1mag'], row['Colour'])), axis = 1).astype(np.float128)
    if(len(MCLR)<1):
        newdrop.append(source)
    else:
        
        MCLR[str(RLF.LRMC)] = MCLR.apply(lambda row: GetLR2(row['Multi_LR'], Getqmc(row[RLF.OptMagA], row['Colour']), Getnmc(row[RLF.OptMagA], row['Colour'])), axis = 1).astype(np.float128)
                
        MCLR.to_csv(str(RLF.LR) %source, columns = ['LofarRDis', str(RLF.PossRA), str(RLF.PossDEC), str(RLF.IDW), str(RLF.IDP), str(RLF.OptMagP), str(RLF.OptMagA), str(RLF.LRMC)], header = True, index = False)


source_list=filter_sourcelist(source_list,newdrop)
        
# For each source in the list find the maximum combined LR and store all the information

def FindMax(source):
    info = pd.read_csv(RLF.LR %source, header = 0, usecols = [str(RLF.PossRA), str(RLF.PossDEC), str(RLF.IDW), str(RLF.IDP), str(RLF.LRMC)])
    info[str(RLF.IDW)] = info[str(RLF.IDW)].map(lambda x: x.strip('b').strip("''"))
    info[str(RLF.IDP)] = info[str(RLF.IDP)].map(lambda x: x.strip('b').strip("''"))
    #info[str(RLF.ID3)] = info[str(RLF.ID3)].map(lambda x: x.strip('b').strip("''"))
    CP = info.loc[info[str(RLF.LRMC)].idxmax()].copy()
    CP['PossFail'] = np.where(CP[str(RLF.LRMC)] < RLC.Lth, 1, 0)
    CP[str(RLF.LSN)] = source
    
    return CP

print('Finding max LRs and writing table')

PossHosts = pd.concat([FindMax(s['Source_Name']) for s in source_list], ignore_index = True, axis = 1)
PossHosts.columns = PossHosts.loc[str(RLF.LSN)]
PossHosts = PossHosts.drop(index = [str(RLF.LSN)])
PossHostsT = PossHosts.transpose()
PossHostsT.to_csv(RLF.PossHosts, header = True, index = True)#,  columns = [str(RLF.LSN), str(RLF.PossRA), str(RLF.PossDEC), str(RLF.IDW), str(RLF.IDP), str(RLF.ID3), str(RLF.OptMagP), str(RLF.OptMagA)]


# In[ ]:


# Number of sources with a 0 max LR and therefore would possibly be a failed LR
# or defined by being closest to LOFAR
print('Total number of fails is',np.sum(PossHostsT['PossFail']))

