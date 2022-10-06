#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 10:56:42 2019

This file contains the functions that allow a user to create a table of all
the known sources in the catalogue.  To reduce this catalogue down to a 
sub-catalogue of nearby sources to the sources which have succuessfully drawn
a ridgeline and create a table.  Then to create text files of the possible
sources in a cutout (size generated and used in the Ridgeline code).

These text files are used to determine the pixel distance of all these possible
sources from the LOFAR catalogue position of the corresponding source, the
perpendicular distance of the possible sources from the corresponding
ridgeline and the distance of the source along the ridgeline. This information 
is stored in a further text file which is used to study the N closest sources 
to the LOFAR catalogue position.

Functions to create tables from catalogues with the required information.
Functions to create .txt files with the catalogue information in.
Functions to create .txt files with the LOFAR distance and Perpendicular
distances in.
Functions to obtain the closest N sources, how many are host ID's and the
distances of all these sources.

This file deals with just the host sources of the ridgelines.  It contains
the functions that are used to determine the perpendicular distance of the
host from the ridgeline of the source and from the LOFAR catalogue position
corresponding to the source.  These are looped so all hosts distances are
contained in two seperate lists.

Functions for obtaining which sources successfully drew ridgelines.
Functions for obtaining perpendicular pixel distances.
Functions for obtaining pixel distance from LOFAR catalogue position.
Functions for obtaining the distance along the ridgeline of the host source.

Edited from 22/01/20

This file contains further functions for generating the perpendicular distance
from the ridgeline and the distance from the LOFAR Catalogue position using
the distance including uncertainties formula in Best et. al. (2003)
(see end of file for ref.), and those functions required to use these distances
to generate the likelihood ratio used in this paper.  All of this information
is stored in text files created in the running of the functions.

Functions to create tables and .txt files.
Functions to calculate R distance to ridgeline and LOFAR Catalogue postion
and create.txt files.
Functions to calculate the corresponding likelihood ratios and create .txt
files.

These functions deal with the host for the radio sources and with the closest 
n sources to the LOFAR catalogue position.

@author: bonnybarkus

"""

#import math
import numpy as np
from numpy.linalg import norm, solve
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS, utils
from ridge_toolkitDR2 import GetProblematicSources
import RLConstantsDR2 as RLC
import RidgelineFilesDR2 as RLF

#############################################

def GetCatArea(lcat):
    allcat=Table.read(lcat)
    ramax=np.nanmax(allcat['RA'])
    ramin=np.nanmin(allcat['RA'])
    decmin=np.nanmin(allcat['DEC'])
    decmax=np.nanmax(allcat['DEC'])
    dd=decmax-decmin
    dmid=(decmin+decmax)/2.0
    dr=(ramax-ramin)*np.cos(dmid*np.pi/180.0)
    a=dr*dd*3600*3600.0
    return a

def AllDistFromLofar(source_list):
    
    """
    Creates a list of all the distances of all possible sources in all the
    cutouts from the LOFAR catalogue position for the corresponding source.
    This list is generated using the .txt files. Once these have been created
    the CreateCutOutCat function does not need to be run again.
    
    Parameters
    ----------
    
    source_list - list,
                  a list of the sources that succesfully drew a ridgeline.    

    Returns
    -------
    
    alldistfromlofar - list,
                       a list of all the distances of all the possible 
                       sources in the cutouts from the LOFAR catalogue
                       position of the corresponding source.
               
    """
    
    alldistfromlofar = []
    
    for source in source_list:
        
        file3 = open(RLF.coc %source, 'r')
        distfromlofar = []
        
        for row in file3:        
            row_info = row.strip()
            row_cols = row_info.split(',')
            dist = float(row_cols[4])
            distfromlofar.append(dist)
        
        for i in range(len(distfromlofar)):
            alldistfromlofar.append(distfromlofar[i])
            
    file3.seek(0) 
    
    return alldistfromlofar

#############################################
    
def CatPerpDisList(source, cutoutcat):
    
    """
    Creates a list of perpendicular distances of all possible sources in 
    the cutout region from the drawn ridgeline for a source.  It does this
    straight from the cutout catalogue so this catalogue needs to be
    generated first.
    
    Parameters
    ----------
    
    source - str,
             the LOFAR ID for the successfully drawn source.
             
    cutoutcat - table,
                a table containing the AllWISE ID, object ID, RA and DEC,
                x and y position in pixels and distance in pixels from the
                LOFAR catalogue position.

    Returns
    -------
    
    pdislist - list,
               an list of the perpendicular distances of all the possible
               sources from the ridgeline of the source.
               
    """
    
    file1 = open(RLF.R1 %source, 'r')
    file2 = open(RLF.R2 %source, 'r')
    
    pdislist = []
    
    points = GetPointList(file1, file2)[0]
    
    for i in range(np.shape(cutoutcat)[0]):
        closest1, closest2 = ClosestPoints(points, cutoutcat['xpix'][i], cutoutcat['ypix'][i])
        m, c = ClosestPointEq(closest1, closest2)
        pdis = PerpDistance(cutoutcat['xpix'][i], cutoutcat['ypix'][i], m, c)
        pdislist.append(pdis)

    return pdislist

#############################################
    
def ClosestPointEq(closest1, closest2):
    
    """
    Returns the gradient and intercept of the equation of the line that
    connects the two points that are closest to the host ID.
    
    Parameters
    ----------
    
    closest1 - array, shape (2, ),
               the closest point to the host ID.
    
    closest2 - array, shape (2, ),
               the second closest point to the host ID.
               
    Returns
    -------
    
    m - float,
        the value of the gradient of the equation of the line joining
        the two closest points.
    
    c - float,
        the value of the intercept of the equation of the line joining
        the two closest points.
        
    """
    
    m, c = np.polyfit([closest1[0], closest2[0]], [closest1[1], closest2[1]], 1)
    
    return np.float(m), np.float(c)

#############################################
    
def ClosestPoints(points, host):
    
    """
    Returns the coordinates of the two closest points on the ridgeline
    to the host ID, and their index in the list of ridgeline points.
    
    Parameters
    ----------
    
    points - list, 
             list of all the points along the ridgeline.
             
    hostx - float,
            the x pixel coordinate for the host ID.
            
    hosty - float,
            the y pixel coordinate for the host ID.
    
    Returns
    -------
    
    closest1 - array, shape (2, ),
               the closest point to the host ID.
               
    closest2 - array, shape (2, ),
               the second closest point to the host ID.
    
    c1index - int,
              the index of the closest point to the host ID.
    
    c2index - int,
              the index of the second closest point to the host ID.
    
    """

    diffs = np.zeros(3)
    
##  Creates the list of differences in the points and the host
    for point in points:
            diff = norm(point - host)
            diff_summary = [float(diff), float(point[0]), float(point[1])]
            diffs = np.vstack((diff_summary, diffs))      
##  Orders that list by smallest difference then finds the smallest and index in the list of points
    diffs = diffs[:-1, :]
    sortdiffs = diffs[diffs[:,0].argsort()]
    closest1 = sortdiffs[0, 1:]
    pointsl  = [point.tolist() for point in points]
    c1index = pointsl.index(list(closest1))
        
##  If at the beginning of the list make the second point the next one in the list
    if c1index == 0:
        c2index = c1index + 1
        closest2 = points[c2index]
##  If at the end of the list make the second one the one before in the list
    elif c1index == len(pointsl) - 1:
        c2index = c1index - 1
        closest2 = points[c2index]
##  Otherwise label each one either side as the next two closest
    else:
        c2index = c1index - 1
        c3index = c1index + 1
        closest2 = points[c2index]
        closest3 = points[c3index]
##  Work out the direction vectors from the closest to the host and to each other point     
        d1 = np.array([host[0] - closest1[0], host[1] - closest1[1]])  ##  Direction to host
        d2 = np.array([closest2[0] - closest1[0], closest2[1] - closest1[1]])  ##  Direction to c2
        d3 = np.array([closest3[0] - closest1[0], closest3[1] - closest1[1]])  ##  Direction to c3
##  Work out the dot product (angle) between the host and each closest
        dot1 = np.dot(d1, d2)  ## Angle between host and closest2
        dot2 = np.dot(d1, d3)  ## Angle between host and closest3
##  Find the one with the acute angle (dot > 0)
        dot = max(dot1, dot2)
##  If this is dot2 mark second point as having the third point information
        if dot == dot2:
            closest2 = points[c3index]
            c2index = c3index
        
    return closest1, closest2, c1index, c2index

#############################################

def CreateAllPerpDisList(source_list):
    
    """
    Creates a list of all the perpendicular distances of all possible sources
    in all the cutouts from the drawn ridgelines for all sources. This list
    is generated using the .txt files. Once these have been created the
    CreateCutOutCat function does not need to be run again.
    
    Parameters
    ----------
    
    source_list - list,
                  a list of the sources that succesfully drew a ridgeline.    

    Returns
    -------
    
    allperpdis - list,
                 a list of all the perpendicular distances of all the 
                 possible sources in the cutouts from the drawn ridgelines 
                 of the sources.
           
    """
    
    allperpdis = []
    
    for source in source_list:
        
        file1 = open(RLF.R1 %source, 'r')
        file2 = open(RLF.R2 %source, 'r')
        file3 = open(RLF.coc %source, 'r')
           
        perpdis = []
            
        points = GetPointList(file1, file2)[0]
            
        for row in file3:
            row_info = row.strip()
            row_cols = row_info.split(',')
            xpix = float(row_cols[2])
            ypix = float(row_cols[3])
            pixloc = np.array([xpix, ypix])
            
            closest1, closest2 = ClosestPoints(points, pixloc)[:2]
            #m, c = ClosestPointEq(closest1, closest2)
            pdis = PerpDistance(closest1, closest2, pixloc)
            perpdis.append(pdis)
        
        file3.seek(0)
        
        for i in range(len(perpdis)):
            allperpdis.append(perpdis[i])
        
    return allperpdis

#############################################

def CreateCutOutCat(source, lofar_table, subcat, Lra, Ldec, lengthpix):
    
    """
    Creates a catalogue of all the possible sources in the region of the
    cutout.  This is within the size of the source from the LOFAR catalogue
    position.  It list all the source's AllWISE ID's, object ID's, RA and DEC
    of the possible source, the x and y position of the possible source in
    pixels, and the distance of the possible source from the LOFAR catalogue
    position in pixels.  This catalogue can be given as an output and is
    saved as a .txt file to '/CutOutCats/' for use with other functions.
    
    Parameters
    ----------
    
    source - str,
             the LOFAR ID for the successfully drawn source.
             
    lofar_table - table,
                  the table containing all the LOFAR data from the .fits
                  catalogue.
                
    subcat - table,
             a table with the same columns as the main .fits catalogue which
             contains only the sources within the given RA and DEC distances
             from the LOFAR catalogue position for the given source.
    
    Constants
    ---------
    
    rdel - float,
           the equivalent pixel value to RA in a FITS file, set value.
           
    ddel - float,
           the equivalent pixel value to DEC in a FITS file, set value.       
        
        
    Returns
    -------
    
    cutoutcat - table,
                a table containing the AllWISE ID, object ID, RA and DEC,
                x and y position in pixels and distance in pixels from the
                LOFAR catalogue position.
             
    """

    #Lra, Ldec, size = SourceInfo(source, lofar_table)[:3]
    
    hdu = fits.open(RLF.fits + source + '.fits')
    header = hdu[0].header
    wcs = WCS(header)
    
    #lengthpix = (size/RLC.ddel/3600.0)  ##  Find the length in a positive value of pixels
    
    sourcecoords = SkyCoord((subcat[str(RLF.PossRA)])*u.degree, (subcat[str(RLF.PossDEC)])*u.degree, \
                            frame = 'fk5')
    sourcepix = utils.skycoord_to_pixel(sourcecoords, wcs, origin = 0)
    
    subcat['xpix'] = sourcepix[0]  ##  Create new columns of the information
    subcat['ypix'] = sourcepix[1]
    
    lofarcoords = SkyCoord(Lra*u.degree, Ldec*u.degree, \
                            frame = 'fk5')
    lofarpix = utils.skycoord_to_pixel(lofarcoords, wcs, origin = 0)
    LofarRa = lofarpix[0]
    LofarDec = lofarpix[1]
    
    ##  Working out everything that is the source size around the lofar position
    cutoutcatx = subcat[subcat['xpix'] < LofarRa + lengthpix]
    cutoutcatx2 = cutoutcatx[cutoutcatx['xpix'] > LofarRa - lengthpix]
    cutoutcaty = cutoutcatx2[cutoutcatx2['ypix'] < LofarDec + lengthpix]
    cutoutcat = cutoutcaty[cutoutcaty['ypix'] > LofarDec - lengthpix]
    
    cutoutcat['xpix'] = cutoutcat['xpix'] - (LofarRa - lengthpix)
    cutoutcat['ypix'] = cutoutcat['ypix'] - (LofarDec - lengthpix)
    
    ##  Calculating the distance from lofar and creating a neew column
    disfromlofar = np.sqrt((cutoutcat['xpix'] - lengthpix) ** 2 + \
                                           (cutoutcat['ypix'] - lengthpix) ** 2)
    
    cutoutcat['disfromLOFAR'] = disfromlofar
    print("Source is "+str(source)+" and length of cat is "+str(len(cutoutcat)))
    lcat=len(cutoutcat)
    np.savetxt(RLF.coc %source, cutoutcat[str(RLF.PossRA), str(RLF.PossDEC), 'xpix', 'ypix', 'disfromLOFAR', 'raErr', 'decErr'], delimiter = ',', fmt='%s', encoding = 'utf-8')

    return cutoutcat,lcat

#############################################
    
def CreateDistFromLofar(source):
    
    """
    Creates a list of distances of all possible sources in the cutout region
    from the LOFAR catalogue postition for the corresponding source.  Creates
    an array of all distances and their AllWISE ID's of all the possible
    sources in the cutout region from the LOFAR catalogue position for the 
    corresponding source.  These are created using the .txt files.  Once 
    these have been created the CreateCutOutCat function does not need to be
    run again.
    
    Parameters
    ----------
    
    source - str,
             the LOFAR ID for the successfully drawn source.
             

    Returns
    -------
    
    distfromlofar - list,
                    a list of the distances of all the possible sources 
                    from the LOFAR catalogue position of the corresponding 
                    source.

    distfromlofarids - array,
                       an array of the distances of the possible sources in
                       the cutout from the LOFAR catalogue position and their
                       corresponding AllWISE ID's.  UNless ID's are needed use
                       distfromlofar, possible accuracy loss due to conversion
                       to string in array
      
    """
    distfromlofar = []
    #distfromlofarids = np.zeros(2)   
   
    file3 = open(RLF.coc %source, 'r')
            
    for row in file3:        
        row_info = row.strip()
        row_cols = row_info.split(',')
        #poss_id = row_cols[0].strip("''").strip('b').strip("''")
        dist = float(row_cols[4])
        distfromlofar.append(dist)
        #distfromlofarids = np.vstack(((poss_id, dist), distfromlofarids))
    
    file3.seek(0)
    
    return distfromlofar, #distfromlofarids[:-1]
    
#############################################
    
def CreateDistTable(source):
    
    """
    Creates a table of distances of all possible sources in a cutout
    catalogue.  The table contains the AllWISE ID, the object ID, the
    distance from the LOFAR catalogue position, the perpendicular distance
    from the ridgeline and distance along the ridgeline (starting from the 
    beginning of ridgeline 1) as a percentage of the total length for the 
    corresponding source.  The output is the table and this is also saved 
    as a .txt file in the directory '/Distances/' for use later.  This 
    function does not need to be re-run, once the .txt files have been 
    created.
    
    Parameters
    ----------
    
    source - str,
             the LOFAR ID for the successfully drawn source.    

    Return
    ------
    
    distances - table,
                a table containing the AllWISE ID, object ID, distance from
                LOFAR catalogue position, perpendicular distance and distance 
                along the ridgeline for every possible source in the cutout 
                for source.
    
    NOTE: To obtain all distance tables for all the sources in the source_list
    this function has to be looped through in the notebook.
    
    """
    
    file1 = open(RLF.R1 %source, 'r')
    file2 = open(RLF.R2 %source, 'r')
    file3 = open(RLF.coc %source, 'r')
    
    points, rllength = GetPointList(file1, file2)
    
    #allwise = []
    #obj_id = []
    lofar_dis = []
    perp_dis = []
    distalong = []
    
    for row in file3:
        row_info = row.strip()
        row_cols = row_info.split(',')
        #allw = row_cols[0].strip("''").strip('b').strip("''")
        #objid = float(row_cols[1])
        xpix = float(row_cols[2])
        ypix = float(row_cols[3])
        ldis = float(row_cols[4])
        pixloc = np.array([xpix, ypix])
        
        closest1, closest2, c1index, c2index = ClosestPoints(points, pixloc)
        pdis = PerpDistance(closest1, closest2, pixloc)
        
        posindex = min(c1index, c2index)
        Ix, Iy = PointOfIntersection(closest1, closest2, pixloc)
        hostlength, totallength = LengthOnLine(rllength, points, posindex, Ix, Iy)
        alongdist = (hostlength/totallength) * 100
        
        
        #allwise.append(allw)
        #obj_id.append(objid)
        lofar_dis.append(ldis)
        perp_dis.append(pdis)
        distalong.append(alongdist)
        content = np.column_stack((lofar_dis, perp_dis, distalong))        
        columns = ['LOFAR_dis', 'Perp_dis', 'Dist Along']
        distances = Table(content, names = columns, dtype = ('f8', 'f8', 'f8'))
        
        np.savetxt(RLF.Dists %source, distances['LOFAR_dis', 'Perp_dis', 'Dist Along'], delimiter = ',', fmt='%s', encoding = 'utf-8')

    print("Written dists file for source "+str(source))
        
    file3.seek(0)

    #return distances
    # may need to add this back, but it fails if coc empty

#############################################
    
def CreatePerpDistList(source):
    
    """
    Creates a list of perpendicular distances of all possible sources in 
    the cutout region from the drawn ridgeline for a source.  Creates an
    array of all perpendicular distances and their AllWISE ID's of
    all the possible sources in the cutout region from the drawn ridgeline
    for a source.  These are created using the .txt files.  Once these have
    been created the CreateCutOutCat function does not need to be run again.
    
    Parameters
    ----------
    
    source - str,
             the LOFAR ID for the successfully drawn source.
             

    Returns
    -------
    
    pdislist - list,
               a list of the perpendicular distances of all the possible
               sources from the ridgeline of the source.

    perpdisids - array,
                 an array of the perpendicular distances of the possible
                 sources in the cutout from their ridgelines and their
                 corresponding AllWISE ID's. Unless ID's are needed use 
                 pdislist, possible loss of accuracy due to string conversion
                 in array.
                 
    """
    
    file1 = open(RLF.R1 %source, 'r')
    file2 = open(RLF.R2 %source, 'r')
    file3 = open(RLF.coc %source, 'r')
        
    perpdis = []
    #perpdisids = np.zeros(2)
        
    points = GetPointList(file1, file2)[0]
        
    for row in file3:
        row_info = row.strip()
        row_cols = row_info.split(',')
        #poss_id = row_cols[0].strip("''").strip('b').strip("''")
        xpix = float(row_cols[2])
        ypix = float(row_cols[3])
        pixloc = np.array([xpix, ypix])
        
        closest1, closest2 = ClosestPoints(points, xpix, ypix)[:2]
        pdis = PerpDistance(closest1, closest2, pixloc)
        perpdis.append(pdis)
        #perpdisids = np.vstack(((poss_id, float(pdis)), perpdisids))
        
    return perpdis, #perpdisids[:-1]

#############################################

def CreatePositionTable(source_list, available_sources, lofar_table):
    
    """
    ***
    ISN'T NEEDED FOR DR2 AND WILL NOT WORK
    ***
    Creates a table listing each soure that successfully draw a ridgeline. 
    The table contains Source name, AllWISE ID, the RA and DEC of the LOFAR 
    Catlogue position, the radio error on the LOFAR Catalogue RA and DEC, 
    the RA and DEC of the host, and the optical errors on the RA and DEC of 
    the host for each source.  The output is the table and this is also 
    saved as a .txt file in the directory '/Distances/' for use later.  This 
    function does not need to be re-run, once the .txt file has been 
    created.   
    
    Parameters
    ----------
    
    source_list - list,
                  a list of the sources that succesfully drew a ridgeline.
                  
    available_sources - 2D array, shape(x, 14), dtype='str'
                        array containing x source names and LOFAR catalogue
                        locations in pixel values, number of components
                        associated with source and LOFAR catalogue locations
                        as RA and DEC, the length of the source in pixels
                        and the total flux for all files available in a 
                        given directory. The host location in pixels, the
                        error on the RA and DEC of the host location, and
                        the RA and DEC of the host.
    
    lofar_table - table,
                  the table containing all the LOFAR data from the .fits
                  catalogue.                       
    
    Returns
    -------
    
    positions - table,
                a table containing the Source name, AllWISE ID, the RA 
                and DEC of the LOFAR Catlogue position, the radio error on
                the LOFAR Catalogue RA and DEC, the RA and DEC of the host,
                and the optical errors on the RA and DEC of the host for each 
                source.
    
    """
    
    sname = []
    AllWISE = []
    lra = []
    ldec = []
    lraerr = []
    ldecerr = []
    hra = []
    hdec = []
    hraerr = []
    hdecerr = []
    
    for source in source_list:
        for asource in available_sources:
            if source == asource[0]:
                source_name = asource[0]
                LRA = float(asource[4])
                LDEC = float(asource[5])
                LRA_errRad = float(asource[10])
                LDEC_errRad = float(asource[11])
                host_RA = float(asource[12])
                host_DEC = float(asource[13])
                
                allwise = SourceInfo(source_name, lofar_table)[3].strip('\n').strip()
                file3 = open(RLF.coc %source_name, 'r')
                                
                for row in file3:
                    row_info = row.strip()
                    row_cols = row_info.split(',')
                    allw = row_cols[0].strip("''").strip('b').strip("''")
                    
                    if allw == allwise:
                        host_RA_errOpt = float(row_cols[7])
                        host_DEC_errOpt = float(row_cols[8])
                
                sname.append(source_name)
                lra.append(LRA)
                ldec.append(LDEC)
                lraerr.append(LRA_errRad)
                ldecerr.append(LDEC_errRad)
                hra.append(host_RA)
                hdec.append(host_DEC)
                hraerr.append(host_RA_errOpt)
                hdecerr.append(host_DEC_errOpt)
                AllWISE.append(allwise)
                
                file3.seek(0)

    content = np.column_stack((sname, AllWISE, lra, ldec, lraerr, ldecerr, hra, hdec, hraerr, hdecerr))        
    columns = ['Source_Name', 'AllWISE', 'LOFAR_RA', 'LOFAR_DEC', 'LOFAR_RA_errRad', 'LOFAR_DEC_errRad', 'Host_RA', 'Host_DEC', 'Host_RA_errOpt', 'Host_DEC_errOpt']
    positions = Table(content, names = columns, dtype = ('S100', 'S100', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))
    posdf = positions.to_pandas() 
    #np.savetxt(RLF.Position, positions['Source_Name', 'AllWISE', 'LOFAR_RA', 'LOFAR_DEC', 'LOFAR_RA_errRad', 'LOFAR_DEC_errRad', 'Host_RA', 'Host_DEC', 'Host_RA_errOpt', 'Host_DEC_errOpt'], delimiter = ',', fmt='%s', encoding = 'utf-8')
    posdf.to_csv(str(RLF.Position), columns = columns, header = True, index = False)
    
    return positions

#############################################

def CreateLDistTable(source, available_sources):
    
    """
    Creates a table of the calculated R distances of all possible sources 
    in a cutout catalogue.  The table contains the AllWISE ID, the calculated 
    sigma for RA and DEC, the R distance from LOFAR catalogue position, the RA
    and DEC, and errors on the RA and DEC for every possible optical source 
    in the cutout for each source.  The output is the table and this is also 
    saved as a .txt file in the directory '/Distances/' for use later.  This 
    function does not need to be re-run, once the .txt files have been 
    created.
    
    Parameters
    ----------
    
    source_list - list,
                  a list of the sources that succesfully drew a ridgeline.

    available_sources - 2D array, shape(x, 14), dtype='str'
                        array containing x source names and LOFAR catalogue
                        locations in pixel values, number of components
                        associated with source and LOFAR catalogue locations
                        as RA and DEC, the length of the source in pixels
                        and the total flux for all files available in a 
                        given directory. The host location in pixels, the
                        error on the RA and DEC of the host location, and
                        the RA and DEC of the host.

    Return
    ------
     
    rdistances - table,
                 a table containing the AllWISE ID, the calculated sigma 
                 for RA and DEC, the R distance from LOFAR catalogue position,
                 the RA and DEC, and errors on the RA and DEC for every 
                 possible optical source in the cutout for each source.
    
    NOTE: To obtain all rdistance tables for all the sources in the source_list
    this function has to be looped through in the notebook.
    
    """
     
    #allwise = []
    lofar_rdis = []
    sigra = []
    sigdec = []
    optRA = []
    optDEC = []
    optRAerr = []
    optDECerr = []
    
    #for source in source_list:
    for asource in available_sources:
        if source == asource['Source_Name']:
            #source_name = asource[0]
            LOFAR_RA = asource['RA']
            LOFAR_DEC = asource['DEC']
            LOFAR_RA_errRad = asource['E_RA']
            LOFAR_DEC_errRad = asource['E_DEC']
            
            file1 = open(RLF.R1 %source, 'r')
            file2 = open(RLF.R2 %source, 'r')
            file3 = open(RLF.coc %source, 'r')
    
            points, rllength = GetPointList(file1, file2)
            
            for row in file3:
                row_info = row.strip()
                row_cols = row_info.split(',')
                #allw = row_cols[0].strip("''").strip('b').strip("''")
                poss_RA = float(row_cols[0])
                poss_DEC = float(row_cols[1])
                poss_RA_errOpt = float(row_cols[5])
                poss_DEC_errOpt = float(row_cols[6])
                
                sigRA, sigDEC = SigmaR(LOFAR_RA_errRad, LOFAR_DEC_errRad, poss_RA_errOpt, poss_DEC_errOpt)
                delRA, delDEC = DeltaRADEC(LOFAR_RA, LOFAR_DEC, poss_RA, poss_DEC)
                r = LDistance(delRA, delDEC)
                                
                #allwise.append(allw)
                sigra.append(sigRA)
                sigdec.append(sigDEC)
                lofar_rdis.append(r)
                optRA.append(poss_RA)
                optDEC.append(poss_DEC)
                optRAerr.append(poss_RA_errOpt)
                optDECerr.append(poss_DEC_errOpt)
                
                content = np.column_stack((sigra, sigdec, lofar_rdis, optRA, optDEC, optRAerr, optDECerr))#, perp_rdis))        
                columns = ['sigRA', 'sigDEC', 'LOFAR_Rdis', 'Opt_RA', 'Opt_DEC', 'Opt_RA_err', 'Opt_DEC_err']#, 'Perp_rdis']
                ldistances = Table(content, names = columns, dtype = ('f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))#, 'f8'))
                
                np.savetxt(RLF.LDists %source, ldistances['sigRA', 'sigDEC', 'LOFAR_Rdis', 'Opt_RA', 'Opt_DEC', 'Opt_RA_err', 'Opt_DEC_err'], delimiter = ',', fmt='%s', encoding = 'utf-8')
            
            file3.seek(0)
       
    #return ldistances

#############################################

def CreateRDistTable(source, available_sources):
    
    """
    Creates a table of the calculated R distances of all possible sources 
    in a cutout catalogue.  The table contains the AllWISE ID, the calculated 
    sigma for RA and DEC, the R distance from LOFAR catalogue position, the RA
    and DEC, and errors on the RA and DEC for every possible optical source 
    in the cutout for each source.  The output is the table and this is also 
    saved as a .txt file in the directory '/Distances/' for use later.  This 
    function does not need to be re-run, once the .txt files have been 
    created.
    
    Parameters
    ----------
    
    source_list - list,
                  a list of the sources that succesfully drew a ridgeline.

    available_sources - astropy table

    Return
    ------
     
    rdistances - table,
                 a table containing the AllWISE ID, the calculated sigma 
                 for RA and DEC, the R distance from LOFAR catalogue position,
                 the RA and DEC, and errors on the RA and DEC for every 
                 possible optical source in the cutout for each source.
    
    NOTE: To obtain all rdistance tables for all the sources in the source_list
    this function has to be looped through in the notebook.
    
    """
     
    #allwise = []
    lofar_rdis = []
    sigra = []
    sigdec = []
    optRA = []
    optDEC = []
    optRAerr = []
    optDECerr = []
    
    #for source in source_list:
    for asource in available_sources:
        if source == asource[0]:
            #source_name = asource[0]
            LOFAR_RA = float(asource[4])
            LOFAR_DEC = float(asource[5])
            LOFAR_RA_errRad = float(asource[8])
            LOFAR_DEC_errRad = float(asource[9])
            
            file1 = open(RLF.R1 %source, 'r')
            file2 = open(RLF.R2 %source, 'r')
            file3 = open(RLF.coc %source, 'r')
    
            points, rllength = GetPointList(file1, file2)
            
            for row in file3:
                row_info = row.strip()
                row_cols = row_info.split(',')
                #allw = row_cols[0].strip("''").strip('b').strip("''")
                poss_RA = float(row_cols[0])
                poss_DEC = float(row_cols[1])
                poss_RA_errOpt = float(row_cols[5])
                poss_DEC_errOpt = float(row_cols[6])
                
                sigRA, sigDEC = SigmaR(LOFAR_RA_errRad, LOFAR_DEC_errRad, poss_RA_errOpt, poss_DEC_errOpt)
                delRA, delDEC = DeltaRADEC(LOFAR_RA, LOFAR_DEC, poss_RA, poss_DEC)
                r = RDistance(delRA, delDEC, sigRA, sigDEC)
                                
                #allwise.append(allw)
                sigra.append(sigRA)
                sigdec.append(sigDEC)
                lofar_rdis.append(r)
                optRA.append(poss_RA)
                optDEC.append(poss_DEC)
                optRAerr.append(poss_RA_errOpt)
                optDECerr.append(poss_DEC_errOpt)
                
                content = np.column_stack((sigra, sigdec, lofar_rdis, optRA, optDEC, optRAerr, optDECerr))#, perp_rdis))        
                columns = ['sigRA', 'sigDEC', 'LOFAR_Rdis', 'Opt_RA', 'Opt_DEC', 'Opt_RA_err', 'Opt_DEC_err']#, 'Perp_rdis']
                rdistances = Table(content, names = columns, dtype = ('f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))#, 'f8'))
                
                np.savetxt(RLF.RDists %source, rdistances['sigRA', 'sigDEC', 'LOFAR_Rdis', 'Opt_RA', 'Opt_DEC', 'Opt_RA_err', 'Opt_DEC_err'], delimiter = ',', fmt='%s', encoding = 'utf-8')
                
            file3.seek(0)
            
    #return rdistances

#############################################

def CreateSubCat(source_table, lofarra, lofardec):
    
    """
    Creates a sub-catalogue from the main catalogue of all the sources
    that lie within a given distance from the given postion.  The intention
    is to use it centred on the LOFAR catalogue postion with set rsize and
    dsize for each source.
    
    Parameters
    ----------
    
    source_table - table,
                   the table containing all the data from the .fits catalogue.
    
    lofarra - float,
              the LOFAR catalogue RA position.
    
    lofardec - float,
               the LOFAR catalogue DEC position.
    
    Constants
    ---------
    
    rsize - float,
            the seperation distance, in arcmins, of the RA which determines
            the size of the sub-catalogue.
    
    dsize - float,
            the seperation distance, in arcmins, of the DEC which determines
            the size of the sub-catalogue.
    
    Returns
    -------
    
    subcat - table,
             a table with the same columns as the main .fits catalogue which
             contains only the sources within the given RA and DEC distances
             from the LOFAR catalogue position for the given source.
              
    """
    subcat = source_table[(np.abs(source_table[RLF.PossRA] - lofarra) * \
                           np.cos(lofardec * np.pi / 180.0) < RLC.rsize) & \
                                (np.abs(source_table[RLF.PossDEC] - lofardec) < RLC.dsize)]
    
    return subcat

#############################################
    
def DeltaRADEC(LOFAR_RA, LOFAR_DEC, Host_RA, Host_DEC):
    
    """
    Calculates the RA and DEC offset between two points.  This is orginally 
    set up to calculate the offset between the host and the LOFAR Catalogue
    position however it gets used throughout the code to calculate the offset
    between two points, hence the parameter names.
    
    Parameters
    ----------
    
    LOFAR_RA - float,
               the value of the RA coordinate of the LOFAR Catalogue position
               in degrees.
             
    LOFAR_DEC - float,
                the value of the DEC coordinates of the LOFAR Catalogue 
                position in degrees.
    
    Host_RA - float,
              the value of the RA coordinate of the host postion in degrees.
    
    Host_DEC - float,
               the value of the DEC coordinate of the host position in degrees.
    
    Returns
    -------
    
    DelRA - float,
            the value of the offset of the RA coordinates between two points
            in arcseconds.
    
    DelDEC - float,
             the value of the offset of the DEC coordinates between two points
             in arcseconds.
    
    """
    
    LOFAR = SkyCoord(LOFAR_RA, LOFAR_DEC, unit = u.deg)  ##  Turn into a pair of 
    Host = SkyCoord(Host_RA, Host_DEC, unit = u.deg)  ##  coordinates to use
    
    OffRA, OffDEC = (LOFAR).spherical_offsets_to(Host)  ##  Work out the correct distance
    ORA = OffRA.deg
    ODEC = OffDEC.deg ##  Remove the coordinate units
    
    DelRA = ORA * 3600  ##  Convert to arcseconds
    DelDEC = ODEC * 3600    
    
    return DelRA, DelDEC
 
#############################################
    
def GetHostLoc(asource):
    
    """
    ***
    NOT NEEDED IN DR2 AND WILL NOT WORK
    ***
    Returns the x and y coordinate in pixels of the host ID, of the LOFAR
    catalogue postion, and the pixel size of the source.  This is retrieved
    from available_sources after it has been run.
    
    Parameters
    ----------
    
    asource - array, shape (10, )
              the line from available sources that corresponds to a
              successfully drawn ridgeline.

    Returns
    -------
    
    hostx - float,
            the x pixel coordinate for the host ID.
            
    hosty - float,
            the y pixel coordinate for the host ID.
    
    LOFARx - float,
             the x pixel coordinate for the LOFAR catalogue.
            
    LOFARy - float,
             the y pixel coordinate for the LOFAR catalogue.    
    
    lmsize - float,
             the pixel size of the source.
             
    """

    lmsize = float(asource[6])
    LOFARx = float(asource[1])
    LOFARy = float(asource[2])
    hostx = float(asource[8]) - (LOFARx - lmsize)  ## Subtract the size of the cutout
    hosty = float(asource[9]) - (LOFARy - lmsize)  ## to bring to the original location
        
    return hostx, hosty, LOFARx, LOFARy, lmsize

#############################################
    
def GetPointList(file1, file2):
    
    """
    Returns a list of coordinate pairs of points along the ridgeline and 
    a list containing the length steps of the line in pixels.  These are 
    obtained from the .txt files generated when the ridgeline is drawn.
    
    Parameters
    ----------
    
    file1 - file path,
            the path to the files produced during the ridgeline drawing
            process. 
            
    file2 - file path,
            the path to the files produced during the ridgeline drawing
            process.
            
    Returns
    -------
    
    ridge_coord - list,
                  a list of the points along the ridgelines.
    
    rllength - list,
               a list of the length steps from the ridgeline drawing files.

    """
    
    ridge_coord = np.zeros(2)
    rllength = []

    for point1 in file1:
        file1_info = point1.strip()
        info1_cols = file1_info.split()
        file1_xcoord = info1_cols[0]   
        file1_ycoord = info1_cols[1]
        sizeR1 = float(info1_cols[3])
        if (file1_xcoord and file1_ycoord) != 'nan': 
            coords1 = [float(file1_xcoord), float(file1_ycoord)]
            ridge_coord = np.vstack((coords1, ridge_coord))
            rllength.append(sizeR1)
    rllength.reverse()
    ridge_coord = ridge_coord[:-1, :]

    R2len = []
    for point2 in file2:
        file2_info = point2.strip()
        info2_cols = file2_info.split()
        file2_xcoord = info2_cols[0]
        file2_ycoord = info2_cols[1]
        sizeR2 = float(info2_cols[3])
        if (file2_xcoord and file2_ycoord) != 'nan': 
            coords2 = [float(file2_xcoord), float(file2_ycoord)]
            R2len.append(sizeR2)
            if coords2 not in ridge_coord:
                ridge_coord = np.vstack((ridge_coord, coords2))
            
    for i in range(1, len(R2len)):
        rllength.append(R2len[i])

    file1.seek(0)
    file2.seek(0)
    
    return ridge_coord, rllength

#############################################

def filter_sourcelist(t, badlist):
    """
    Removes sources with Source_Name in badlist from a table and returns the table

    Parameters
    ----------

    t - astropy table
    badlist -- a list of source names

    Returns
    -------

    t[filt] - a new table with only the sources not on the bad source list
    """

    filt=[]
    for source in t:
        filt.append(source['Source_Name'] not in badlist)
    return t[filt]

def GetSourceList(available_sources, probfile):

    """
    Returns a list of source names for the sources that drew successful
    ridgelines from the selection criteria.  It uses the functions
    GetAvailableSources and GetProblematicSources from Ridge_toolkitNOID, 
    therefore these need to be available to be accessed.

    Parameters
    ----------

    available_sources - astropy table

    probfile - file path,
               this is the problematic source list that is generated
               by running the ridgeline program on the sources obtained 
               from TotalFluxSelector list.

    Returns
    -------

    filtered_sources: a new table with only the non-problematic sources

    NOTE: Ridge_toolkitNOID has to be in the same folder and the two
    files have to be from the sample that is being graphed. Otherwise this
    will not work and/or the wrong list of sources will be obtained.  Also
    be aware of the choice to use full or cutout FITS files, NOID is 
    for cutouts.

    """
    
    problematic_sources = GetProblematicSources(probfile)

    filt=[]
    for source in available_sources:
        filt.append(source['Source_Name'] not in problematic_sources)

    return filter_sourcelist(available_sources,problematic_sources)

#############################################

def HostDistAlongLine(source_list, available_sources):
    
    """
    ***
    NOT NEEDED IN DR2 AND WILL NOT WORK
    ***
    Returns all the host ID distances along their respective ridgelines.
    It returns a list of just the distances and an array of the distances
    and their respective LOFAR source ID's.  The distances are given as a 
    percentage of the total length of the corresponding ridgeline and are 
    taken from the same starting point of the beginning of ridgeline 1.
    
    Parameters
    ----------
    
    source_list - list,
                   a list of the sources that succesfully drew a ridgeline.
                   
    available_sources - 2D array, shape(x, 14), dtype='str'
                        array containing x source names and LOFAR catalogue
                        locations in pixel values, number of components
                        associated with source and LOFAR catalogue locations
                        as RA and DEC, the length of the source in pixels
                        and the total flux for all files available in a 
                        given directory. The host location in pixels, the
                        error on the RA and DEC of the host location, and
                        the RA and DEC of the host.    
    
    Returns
    -------
    
    distalong - list,
                a list of the  distances of the host ID's along their 
                ridgelines as a percentage.
              
    distalongids - array,
                   an array of the  distances of the host ID's along their 
                   ridgelines as a percentage and their corresponding LOFAR
                   source ID's.  Unless ID's are required distalong should
                   be used in case of loss of precision converting to string 
                   in array.
    
    """
    
    distalong = []
    distalongids = np.zeros(2)
    
    for source in source_list:
        for asource in available_sources:
            if source == asource[0]:
                source_name = asource[0]
                
                file1 = open(RLF.R1 %source_name, 'r')
                file2 = open(RLF.R2 %source_name, 'r')            
                
                points, rllength = GetPointList(file1, file2)
                hostx, hosty = GetHostLoc(asource)[:2]
                host = np.array([hostx, hosty])
                c1, c2, c1index, c2index = ClosestPoints(points, host)
                posindex = min(c1index, c2index)
                  
                Ix, Iy = PointOfIntersection(c1, c2, host)
                
                hostlength, totallength = LengthOnLine(rllength, points, posindex, Ix, Iy)
                
                alongdist = (hostlength/totallength) * 100
                
                distalong.append(alongdist)
                distalongids = np.vstack(((asource[0], alongdist), distalongids))

    return distalong, distalongids[:-1]

#############################################

def HostDisFromLofar(source_list, available_sources):
    
    """
    ***
    NOT NEEDED IN DR2 AND WILL NOT WORK
    ***
    Returns an array of distances of the host IDs from their respective
    LOFAR catalogue positions and an array of distances and the LOFAR
    source ids.
    
    Parameters
    ----------
    
    source_list - list,
                   a list of the sources that succesfully drew a ridgeline.
                   
    available_sources - 2D array, shape(x, 14), dtype='str'
                        array containing x source names and LOFAR catalogue
                        locations in pixel values, number of components
                        associated with source and LOFAR catalogue locations
                        as RA and DEC, the length of the source in pixels
                        and the total flux for all files available in a 
                        given directory. The host location in pixels, the
                        error on the RA and DEC of the host location, and
                        the RA and DEC of the host.
                        
    Returns
    -------
    
    distfromlofar - list,
                    a list of the distances of the host ID's from
                    their LOFAR catalogue positions.
              
    distfromlofarids - array,
                       an array of the distances of the host ID's
                       from their LOFAR catalogue postions and their
                       corresponding LOFAR source ID's.  Unless ID's is 
                       required distfromlofar should be used in case of loss
                       of precision in conversion to string in array.
                 
    """
    
    distfromlofar = []
    distfromlofarids = np.zeros(2)
    
    for source in source_list:
        for asource in available_sources:
            if source == asource[0]:    
                LOFARx = float(asource[1])
                LOFARy = float(asource[2])
                hostx = float(asource[8])
                hosty = float(asource[9])
                dist = np.sqrt((hostx - LOFARx) ** 2 + (hosty - LOFARy) ** 2)
                distfromlofar.append(dist)
                distfromlofarids = np.vstack(((asource[0], dist), distfromlofarids))
    
    return distfromlofar, distfromlofarids[:-1]

#############################################

def HostLikelihoodRatio(table):
    
    """
    ***
    NOT NEEDED IN DR2 AND WILL NOT WORK (AS TABLES DO NOT WORK)
    ***
    Calculates Best et. al.'s likelihood ratio using the R distances from 
    the ridgeline.
    NOTE: Can be used on the appropriate table to calculate the likelihood
    ratios for just the hosts.
    
    Parameters
    ----------
    
    table - array, shape(3,)
            the table of R distances and associated uncertainties on RA 
            and DEC.
    
    Returns
    -------
    
    likelihoodratio - list,
                      a list of likelihood ratios for the R distances from
                      the ridgeline.
    
    LRatDistance - array, shape(2,)
                   an array of the likelihood ratios and the corresponding
                   R distance from the ridgeline.    
    
    """
    
    likelihoodratio = []
    distance = []
    
    for row in table:
    
        r = row[0]
        sigmara = row[1]
        sigmadec = row[2]
        Y = Lambda(sigmara, sigmadec)
        
        LR = (np.float128(1.0) / (np.float128(2.0) * Y)) * np.exp(((r ** np.float128(2.0)) / np.float128(2.0)) * (np.float128(2.0) * Y - np.float128(1.0)))
       
        likelihoodratio.append(LR)
        distance.append(r)
    
    LRatDistance = np.column_stack((distance, likelihoodratio))
        
    return likelihoodratio, LRatDistance

#############################################
    
def HostPerpDis(source_list, available_sources):
    
    """
    ***
    NOT NEEDED IN DR2 AND WILL NOT WORK
    ***
    Returns a list of perpedicular distances of the host IDs from
    their respective ridgelines and an array of perpendicular distances
    and the LOFAR source ids.
    
    Parameters
    ----------
    
    source_list - list,
                   a list of the sources that succesfully drew a ridgeline.
                   
    available_sources - 2D array, shape(x, 14), dtype='str'
                        array containing x source names and LOFAR catalogue
                        locations in pixel values, number of components
                        associated with source and LOFAR catalogue locations
                        as RA and DEC, the length of the source in pixels
                        and the total flux for all files available in a 
                        given directory. The host location in pixels, the
                        error on the RA and DEC of the host location, and
                        the RA and DEC of the host.
                        
    Returns
    -------
    
    perpdis - list,
              a list of the perpendicular distances of the host ID's from
              their ridgelines.
              
    perpdisids - array,
                 an array of the perpendicular distances of the host ID's
                 from their ridgelines and their corresponding LOFAR
                 source ID's. Unless ID's and distance is needed perpdis
                 should be used in case of loss of precision due to string
                 conversion in array.
              
    """
    
    perpdis = []
    perpdisids = np.zeros(2)
    
    for source in source_list:
        for asource in available_sources:
            if source == asource[0]:
                source_name = asource[0]
        
                file1 = open(RLF.R1 %source_name, 'r')
                file2 = open(RLF.R2 %source_name, 'r') 
        
                points = GetPointList(file1, file2)[0]
                hostx, hosty = GetHostLoc(asource)[:2]
                host = np.array([hostx, hosty])
                closest1, closest2 = ClosestPoints(points, host)[:2]
                #m, c = ClosestPointEq(closest1, closest2)
                pdis = PerpDistance(closest1, closest2, host)
                perpdis.append(pdis)
                perpdisids = np.vstack(((asource[0], pdis), perpdisids))
    
    return perpdis, perpdisids[:-1]

#############################################

def HostRFromLofar(positions):
    
    """
    ***
    NOT NEEDED IN DR2 AND WILL NOT WORK
    ***
    Finds the R distance as defined in Best et. al. (2003).  This is the
    distance between two points including the calculated uncertainities on
    the RA and DEC.  The R distance is calculated between the host and the 
    LOFAR Catalogue position for each source of a successfully drawn
    ridgeline.
    
    Parameters
    ----------
    
    positions - table,
                a table containing the Source name, AllWISE ID, the RA 
                and DEC of the LOFAR Catlogue position, the radio error on
                the LOFAR Catalogue RA and DEC, the RA and DEC of the host,
                and the optical errors on the RA and DEC of the host for each 
                source.
    
    Returns
    -------
    
    HostR - list,
            a list of the R distances between the hosts and the LOFAR
            Catalogue position.
    
    HostRrandsig - array, shape( , 3)
                   an array contain the R distances between the host and the
                   LOFAR catalogue position, and the corresponding calculated
                   uncertainities for the RA and DEC of each successfull 
                   source.
    
    """

    HostR = []
    sigmara = []
    sigmadec = []
    
    for row in positions:
        LOFAR_RA = row[2]
        LOFAR_DEC = row[3]
        LOFAR_RA_errRad = row[4]
        LOFAR_DEC_errRad = row[5]
        Host_RA = row[6]
        Host_DEC = row[7]
        Host_RA_errOpt = row[8]
        Host_DEC_errOpt = row[9]
        
        SigRA, SigDEC = SigmaR(LOFAR_RA_errRad, LOFAR_DEC_errRad, Host_RA_errOpt, Host_DEC_errOpt)
        
        DelRA, DelDEC = DeltaRADEC(LOFAR_RA, LOFAR_DEC, Host_RA, Host_DEC)
        
        r = RDistance(DelRA, DelDEC, SigRA, SigDEC)
        
        HostR.append(r)
        sigmara.append(SigRA)
        sigmadec.append(SigDEC)
        HostRrandsig = np.column_stack((HostR, sigmara, sigmadec))
        
    return HostR, HostRrandsig

#############################################

def HostRFromRL(source_list, available_sources, positions):
    
    """
    ***
    NOT NEEDED IN DR2 AND WILL NOT WORK
    ***
    Calculates Best et. al. (2003) R distance of the host from the succesfully
    drawn ridgeline.
    
    Parameters
    ----------
 
    source_list - list,
                  a list of the sources that succesfully drew a ridgeline.         

    available_sources - 2D array, shape(x, 14), dtype='str'
                        array containing x source names and LOFAR catalogue
                        locations in pixel values, number of components
                        associated with source and LOFAR catalogue locations
                        as RA and DEC, the length of the source in pixels
                        and the total flux for all files available in a 
                        given directory. The host location in pixels, the
                        error on the RA and DEC of the host location, and
                        the RA and DEC of the host.    
    
    positions - table,
                a table containing the Source name, AllWISE ID, the RA 
                and DEC of the LOFAR Catlogue position, the radio error on
                the LOFAR Catalogue RA and DEC, the RA and DEC of the host,
                and the optical errors on the RA and DEC of the host for each 
                source.    
    
    Returns
    -------
    
    HostRRL - list,
                a list of the R distances between the hosts and the ridgeline.
    
    RLrandsig - array, shape( , 3)
                   an array contain the R distances between the host and the
                   ridgeline, and the corresponding calculated uncertainities
                   for the RA and DEC of each successfull source.    
    
    """
    
    HostRRL = []
    SigmaRA = []
    SigmaDEC = []
    
    for source in source_list:
        for asource in available_sources:
            if source == asource[0]:
                source_name = asource[0]
                lofarx = float(asource[1])
                lofary = float(asource[2])
                lmsize = float(asource[6])
                hostx = float(asource[8]) - (lofarx - lmsize)
                hosty = float(asource[9]) - (lofary - lmsize)
                hostRA = float(asource[12])
                hostDEC = float(asource[13])
                host = np.array([hostx, hosty])
                
                file1 = open(RLF.R1 %source_name, 'r')
                file2 = open(RLF.R2 %source_name, 'r')
                hdu = fits.open(RLF.fits + source_name + '.fits')
                header = hdu[0].header
                wcs = WCS(header)
                
                points = GetPointList(file1, file2)[0]  ## Retrieve the list of points along the RL
                c1, c2 = ClosestPoints(points, host)[:2]  ##  Find the closest two RL points to the host in pixels
                Ix, Iy = PointOfIntersection(c1, c2, host)  ##  Find the point of intersection along the RL
                
                Ixa = Ix + (lofarx - lmsize)  ##  Adjust back to the main picture
                Iya = Iy + (lofary - lmsize)                
                
                ##  Turn the adjusted point of intersection into degrees instead of pixels
                Ideg = utils.pixel_to_skycoord(float(Ixa), float(Iya), wcs, origin = 0)
                Ira = Ideg.ra.degree
                Idec = Ideg.dec.degree
                
                
                for row in positions:                  
                    if source_name == row[0]:  ## find the optical errors from position_info
                        optRAerr = float(row[8])
                        optDECerr = float(row[9])
                
                radRAerr = np.float128(RLC.radraerr)  ##  define the radio errors as zero on the RL
                radDECerr = np.float128(RLC.raddecerr)  ##  until I know better
                
                RLsigRA, RLsigDEC = SigmaR(radRAerr, radDECerr, optRAerr, optDECerr)
                RLdelRA, RLdelDEC = DeltaRADEC(Ira, Idec, hostRA, hostDEC)
                
                RLr = RDistance(RLdelRA, RLdelDEC, RLsigRA, RLsigDEC)
                
                HostRRL.append(RLr)
                SigmaRA.append(RLsigRA)
                SigmaDEC.append(RLsigDEC)
                
                RLrandsigma = np.column_stack((HostRRL, SigmaRA, SigmaDEC))
                
                file1.seek(0)
                file2.seek(0)
    
    return HostRRL, RLrandsigma

#############################################
    
def HostUnderN(source_list, lofar_table, n):
    
    """
    ***
    ISN'T NEEDED FOR DR2 AND WILL NOT WORK
    ***
    Counts the number of hosts which occur within the Nth closest possible
    sources to the LOFAR catalogue position.  Returns the number of hosts
    and a list of all the distances of the Nth closest sources to the LOFAR 
    catalogue position.
    
    Parameters
    ----------
    
    source_list - list,
                  a list of the sources that succesfully drew a ridgeline.     
    
    lofar_table - table,
                  the table containing all the LOFAR data from the .fits
                  catalogue.    
    
    n - int,
        the number of possible sources away from the LOFAR catalogue position,
        or the Nth closest source to the LOFAR catalogue position.
    
    Returns
    -------
    
    host_counter - float,
                   the number of hosts which are present within the Nth 
                   closest sources to the LOFAR catalogue position.
                 
    alldistundern - list,
                    a list of all the distance from the LOFAR catalogue
                    position of all the sources which are the Nth closest to
                    to the LOFAR catalogue position.
    
    NOTES: This function is void due to NClosestSourcesDistances function
    which finds all the different distances.
    
    """
    
    host_counter = 0
    alldistundern = []

    for source in source_list:
        Ldist = []
        Lids = []
        allwise = SourceInfo(source, lofar_table)[3]  ##  Selects only the allwise ID

        distfromlofarids = CreateDistFromLofar(source)[1]
        
        for i in range(np.shape(distfromlofarids)[0]):
            Ldist.append(float(distfromlofarids[i][1]))
            Lids.append(distfromlofarids[i][0])
        
        Lsortids = [Lids for Ldist, Lids in sorted(zip(Ldist, Lids))]# sorts the ids according to the distances
        Lnsortids = Lsortids[:n]
        Lnsortdist = sorted(Ldist)[:n]
        
        # Check if the host id is in the list
        for j in range(len(Lnsortids)):
            if allwise == Lnsortids[j]:
                host_counter += 1 # Increase the host counter by 1 if it is
                
        for k in range(len(Lnsortdist)):
            alldistundern.append(Lnsortdist[k])
                
    return host_counter, alldistundern  

#############################################

def Lambda(SigRA, SigDEC):
    
    """
    Calculates the variable Lambda from the likelihood ratio formula in 
    Best et. al. (2003).
    
    Parameters
    ----------
    
    SigRA - float,
            the uncertainty on the RA in arcseconds.
            
    SigDEC - float,
             the uncertainty on the DEC in arcseconds.
             
    Constants
    ---------
    
    optcount - int,
               the number of optical sources in the optical catalogue found 
               in the same sky area as the radio catalogue.
    
    catarea - int,
              the area of the sky covered by the catalogue in arcseconds
              squared.
    
    Returns
    -------
    
    lamb - float128,
           the variable lambda to be used in the likelihood ratio calculation.
    
    """
    
    optdensity = RLC.optcount / RLC.catarea  ##  The density of optical sources per arcsec sq
    
    lamb = np.pi * SigRA * SigDEC * optdensity
    
    return np.float128(lamb)

#############################################
    
def LengthOnLine(rllength, points, posindex, Ix, Iy):
    
    """
    Calculates the lengths along a ridgeline.  It works out the length of the
    host along the ridgeline from the start of the drawn ridgeline 1, to 
    where it crosses at a perpendicular angle and it calculates the total 
    length of the ridgeline.
    
    Parameters
    ----------

    rllength - list,
               a list of the length steps from the ridgeline drawing files.
    
    points - list, 
             list of all the points along the ridgeline.
    
    posindex - float,
               the minimum of the two index values of the closest points 
               along the ridgeline to the host ID.
    
    Ix - float,
         the x pixel coordinate of the intercept.
         
    Iy - float,
         the y pixel coordinate of the intercept.
    
    Returns
    -------
    
    hostlength - float,
                 the length, in pixels, along the ridgeline where the host ID
                 meets at a perpendicular angle.
    
    totallength - float,
                  the total length, in pixels, of the ridgeline.
    
    """
    
    baselength = 0
    for i in range(posindex - 1):
        diff = np.abs(rllength[i] - rllength[i + 1])
        baselength += diff 
    
    nearest = points[posindex]
    adist = np.sqrt((nearest[0] - Ix) ** 2 + (nearest[1] - Iy) ** 2)
    hostlength = baselength + adist 
    totallength = rllength[0] + rllength[-1]

    return hostlength, totallength

#############################################

def LikelihoodRatios(available_sources,debug=False):
    
    """
    Calculates the likelihood ratio of the n closest sources as determined by
    the function NClosestRDistances and found in the corresponding .txt files
    for the R distances to the LOFAR Catalogue position and to the ridgeline.
    Two seperate .txt files are created for each source containing the AllWise
    ID, the LOFAR or Ridge R distance, the LOFAR or Ridge Likelihood Ratio and
    coordinates of the source corresponding to these values.
    
    Parameters
    ----------
    
    available_sources - astropy table with sources with successfully drawn ridge lines
    """

    # MJH amended this to work entirely in co-ords of the cutout image
    # -- which should always be centred on the LOFAR RA and Dec --
    # rather than the original FITS image.
    
    for asource in available_sources:
        source_name=asource['Source_Name']
        lofarra = asource['RA']
        lofardec = asource['DEC']
        lmsize = asource['Size']/(3600.0*RLC.ddel) # pixels

        if debug: print('Doing source',source_name,'with lmsize',lmsize)
                
        file1 = open(RLF.R1 %source_name, 'r')
        file2 = open(RLF.R2 %source_name, 'r')
        file7 = open(RLF.NDist %source_name)            
        hdu = fits.open(RLF.fitscutout + source_name + '-cutout.fits') # refer to cutout file should mean pixels are right
        header = hdu[0].header
        if debug: print('cutout image shape is',hdu[0].data.shape)
        wcs = WCS(header)

        
        AllLofarLR = []
        AllRidgeLR = []
        #AllWISE = []
        LofarDist = []
        RidgeDist = []
        PossRA = []
        PossDEC = []

        for row in file7:
            row_info = row.strip()
            row_cols = row_info.split(',')
            #AWISE = row_cols[0].strip("''").strip('b').strip("''")
            LofarRDist = float(row_cols[0])  ## R Distance to LOFAR from file
            #LofarSigRA = float(row_cols[2])  ## Sigma RA from file
            #LofarSigDEC = float(row_cols[3])  ## Sigma DEC from file
            Poss_RA = float(row_cols[3])  ## RA of point for the R dist and sigmas
            Poss_DEC = float(row_cols[4])  ## DEC of point for the R dist and sigmas
            Poss_RA_err = float(row_cols[5])  ## opt err on RA for the point
            Poss_DEC_err = float(row_cols[6])  ## opt err on the DEC for the point
            Poss_Coords = SkyCoord(Poss_RA*u.degree, Poss_DEC*u.degree, \
                        frame = 'fk5')  ##  turn RA and DEC in to a single, degree point
            Poss_pix = utils.skycoord_to_pixel(Poss_Coords, wcs, origin = 0)  ## convert to pixels
            Poss_pixx = Poss_pix[0]
            Poss_pixy = Poss_pix[1]
            if debug: print('Position of optical source in pix is',Poss_pixx,Poss_pixy)
            Opt_Poss = np.array([Poss_pixx, Poss_pixy]) ## Combine to an array

            points = GetPointList(file1, file2)[0]  ## Retrieve the list of points along the RL
            c1, c2 = ClosestPoints(points, Opt_Poss)[:2]  ##  Find the closest two RL points to the host in pixels
            if debug: print('c1 and c2 are',c1,c2)
            Ix, Iy = PointOfIntersection(c1, c2, Opt_Poss)  ##  Find the point of intersection along the RL
            if debug: print('Ix,Iy are',Ix,Iy)

            ##  Turn the adjusted point of intersection into degrees instead of pixels
            Ideg = utils.pixel_to_skycoord(float(Ix), float(Iy), wcs, origin = 0)
            Ira = Ideg.ra.degree
            Idec = Ideg.dec.degree

            if debug:
                print('Poss_RA is',Poss_RA,'and Poss_DEC is',Poss_DEC)
                print('Ira is',Ira)
                print('Idec is',Idec)
            
            radRAerr = np.float128(RLC.radraerr)  ##  define the radio errors as zero on the RL
            radDECerr = np.float128(RLC.raddecerr)  ##  until I know better

            RLsigRA, RLsigDEC = SigmaR(radRAerr, radDECerr, Poss_RA_err, Poss_DEC_err)
            #RLsigRA = np.float128(0.3)
            RLdelRA, RLdelDEC = DeltaRADEC(Ira, Idec, Poss_RA, Poss_DEC)
            if debug: print('RLdelRA is',RLdelRA,'and RLdelDEC is',RLdelDEC)
            RidgeRDist = LDistance(RLdelRA, RLdelDEC)
            if debug: print('RidgeRDist is',RidgeRDist,'arcsec')
            RidgeY = Lambda(RLsigRA, RLsigDEC)

            #optdensity = RLC.optcount / area
            optdensity = 1.0 # test

            #RidgeLR = (optdensity / (np.float128(2.0) * RidgeY)) * np.exp(((RidgeRDist ** np.float128(2.0)) / np.float128(2.0)) * ((np.float128(2.0) * RidgeY) - np.float128(1.0))) # Best et al
            RidgeLR = (np.float128(1.0) / (np.sqrt(2.0*np.pi*(1/optdensity)) * RLsigRA)) * np.exp(((-RidgeRDist ** np.float128(2.0)) /(2.0*RLsigRA*RLsigRA)))
            #RidgeLR = (np.sqrt(optdensity) / (np.sqrt(2.0 * np.pi) * RLsigRA)) * np.exp(((-RidgeRDist ** np.float128(2.0)) /(2.0 * RLsigRA ** np.float128(2.0)))) # 1D Gaussian
            #RidgeLR = (np.float128(1.0) / (np.float128(2.0) * np.pi * (np.float128(RLsigRA) ** np.float128(2.0)))) * np.exp((-RidgeRDist ** np.float128(2.0)) / (np.float(2.0) * (np.float128(RLsigRA) ** np.float128(2.0)))) # 2D Gaussian

            #LofarY = Lambda(LofarSigRA, LofarSigDEC)

            #LofarRDist = LofarRDist * lmsize
            LofarLR = (np.float128(1.0) / (np.float128(2.0) * np.pi * (np.float128(RLC.SigLC) ** np.float128(2.0)))) * np.exp((-LofarRDist ** np.float128(2.0)) / (np.float128(2.0) * (np.float128(RLC.SigLC) ** np.float128(2.0)))) # 2D Gaussian

            #AllWISE.append(AWISE)
            AllLofarLR.append(LofarLR)
            AllRidgeLR.append(RidgeLR)
            LofarDist.append(LofarRDist)
            RidgeDist.append(RidgeRDist)
            PossRA.append(Poss_RA)
            PossDEC.append(Poss_DEC)
                    
        RidgeData = np.column_stack((RidgeDist, AllRidgeLR, PossRA, PossDEC))
        LofarData = np.column_stack((LofarDist, AllLofarLR, PossRA, PossDEC))
                
        Lofarcolumns = ['Lofar_R_Distance', 'Lofar_LR', 'PossRA', 'PossDEC']
        LofarInfo = Table(LofarData, names = Lofarcolumns, dtype = ('f8', 'f8', 'f8', 'f8'))
        Ridgecolumns = ['Ridge_R_Distance', 'Ridge_LR', 'PossRA', 'PossDEC']     
        RidgeInfo = Table(RidgeData, names = Ridgecolumns, dtype = ('f8', 'f8', 'f8', 'f8'))
        LofarInfodf = LofarInfo.to_pandas()
        RidgeInfodf = RidgeInfo.to_pandas()
        
        LofarInfodf.to_csv(RLF.LLR %source_name, columns = Lofarcolumns, header = True, index = False)
        RidgeInfodf.to_csv(RLF.RLR %source_name, columns = Ridgecolumns, header = True, index = False)
        #np.savetxt('Ratios/NearestLofarLikelihoodRatios-%s.txt' %source_name, LofarInfo['AllWISE','Lofar_R_Distance', 'Lofar_LR', 'PossRA', 'PossDEC'], delimiter = ',', fmt='%s', encoding = 'utf-8')    
        #np.savetxt('Ratios/NearestRidgeLikelihoodRatios-%s.txt' %source_name, RidgeInfo['AllWISE','Ridge_R_Distance', 'Ridge_LR', 'PossRA', 'PossDEC'], delimiter = ',', fmt='%s', encoding = 'utf-8')   

#############################################

def SimplifiedLR(source_list, available_sources):
    
    """
    Calculates the likelihood ratio of the n closest sources as determined by
    the function NClosestRDistances and found in the corresponding .txt files
    for the R distances to the LOFAR Catalogue position and to the ridgeline.
    Two seperate .txt files are created for each source containing the AllWise
    ID, the LOFAR or Ridge R distance, the LOFAR or Ridge Likelihood Ratio and
    coordinates of the source corresponding to these values.
    
    Parameters
    ----------
    
    source_list - list,
                  a list of the sources that succesfully drew a ridgeline.         

    available_sources - 2D array, shape(x, 14), dtype='str'
                        array containing x source names and LOFAR catalogue
                        locations in pixel values, number of components
                        associated with source and LOFAR catalogue locations
                        as RA and DEC, the length of the source in pixels
                        and the total flux for all files available in a 
                        given directory. The host location in pixels, the
                        error on the RA and DEC of the host location, and
                        the RA and DEC of the host.      
        
    """

    for source in source_list:
        for asource in available_sources:
            if source == asource[0]:
                source_name = asource[0]
                lofarx = float(asource[1])
                lofary = float(asource[2])
                lmsize = float(asource[6])
                
                file1 = open(RLF.R1 %source_name, 'r')
                file2 = open(RLF.R2 %source_name, 'r')
                file7 = open(RLF.NDist %source_name)            
                hdu = fits.open(RLF.fits + source_name + '.fits')
                header = hdu[0].header
                wcs = WCS(header)
                
                AllLofarLR = []
                AllRidgeLR = []
                #AllWISE = []
                LofarDist = []
                RidgeDist = []
                PossRA = []
                PossDEC = []
                
                for row in file7:
                    row_info = row.strip()
                    row_cols = row_info.split(',')
                    #AWISE = row_cols[0].strip("''").strip('b').strip("''")
                    LofarRDist = float(row_cols[1])  ## R Distance to LOFAR from file
                    LofarSigRA = float(row_cols[2])  ## Sigma RA from file
                    LofarSigDEC = float(row_cols[3])  ## Sigma DEC from file
                    Poss_RA = float(row_cols[4])  ## RA of point for the R dist and sigmas
                    Poss_DEC = float(row_cols[5])  ## DEC of point for the R dist and sigmas
                    Poss_RA_err = float(row_cols[6])  ## opt err on RA for the point
                    Poss_DEC_err = float(row_cols[7])  ## opt err on the DEC for the point
                    Poss_Coords = SkyCoord(Poss_RA*u.degree, Poss_DEC*u.degree, \
                                frame = 'fk5')  ##  turn RA and DEC in to a single, degree point
                    Poss_pix = utils.skycoord_to_pixel(Poss_Coords, wcs, origin = 0)  ## convert to pixels
                    Poss_pixx = Poss_pix[0] - (lofarx - lmsize) ##  Adjust to cutout
                    Poss_pixy = Poss_pix[1] - (lofary - lmsize) ##  Adjust to cutout
                    Opt_Poss = np.array([Poss_pixx, Poss_pixy]) ## Combine to an array
                    
                    points = GetPointList(file1, file2)[0]  ## Retrieve the list of points along the RL
                    c1, c2 = ClosestPoints(points, Opt_Poss)[:2]  ##  Find the closest two RL points to the host in pixels
                    Ix, Iy = PointOfIntersection(c1, c2, Opt_Poss)  ##  Find the point of intersection along the RL
                    
                    Ixa = Ix + (lofarx - lmsize)  ##  Adjust back to the main picture
                    Iya = Iy + (lofary - lmsize)                
                    
                    ##  Turn the adjusted point of intersection into degrees instead of pixels
                    Ideg = utils.pixel_to_skycoord(float(Ixa), float(Iya), wcs, origin = 0)
                    Ira = Ideg.ra.degree
                    Idec = Ideg.dec.degree
                    
                    radRAerr = np.float128(RLC.radraerr)  ##  define the radio errors as zero on the RL
                    radDECerr = np.float128(RLC.raddecerr)  ##  until I know better
                    
                    RLsigRA, RLsigDEC = SigmaR(radRAerr, radDECerr, Poss_RA_err, Poss_DEC_err)
                    RLdelRA, RLdelDEC = DeltaRADEC(Ira, Idec, Poss_RA, Poss_DEC)
                    
                    RidgeRDist = RDistance(RLdelRA, RLdelDEC, RLsigRA, RLsigDEC)
                                     
                    RidgeY = Lambda(RLsigRA, RLsigDEC)
            
                    RidgeLR = (np.float128(1.0) / (np.float128(2.0) * RidgeY)) * np.exp(((RidgeRDist ** np.float128(2.0)) / np.float128(2.0)) * (np.float128(2.0) * RidgeY - np.float128(1.0)))
           
                    LofarY = Lambda(LofarSigRA, LofarSigDEC)
                
                    LofarLR = (np.float128(1.0) / (np.float128(2.0) * LofarY)) * np.exp(((LofarRDist ** np.float128(2.0)) / np.float128(2.0)) * (np.float128(2.0) * LofarY - np.float128(1.0)))
                    
                    #AllWISE.append(AWISE)
                    AllLofarLR.append(LofarLR)
                    AllRidgeLR.append(RidgeLR)
                    LofarDist.append(LofarRDist)
                    RidgeDist.append(RidgeRDist)
                    PossRA.append(Poss_RA)
                    PossDEC.append(Poss_DEC)
                    
        RidgeData = np.column_stack((RidgeDist, AllRidgeLR, PossRA, PossDEC))
        LofarData = np.column_stack((LofarDist, AllLofarLR, PossRA, PossDEC))
                
        Lofarcolumns = ['Lofar_R_Distance', 'Lofar_LR', 'PossRA', 'PossDEC']
        LofarInfo = Table(LofarData, names = Lofarcolumns, dtype = ('f8', 'f8', 'f8', 'f8'))
        Ridgecolumns = ['Ridge_R_Distance', 'Ridge_LR', 'PossRA', 'PossDEC']     
        RidgeInfo = Table(RidgeData, names = Ridgecolumns, dtype = ('f8', 'f8', 'f8', 'f8'))
        LofarInfodf = LofarInfo.to_pandas()
        RidgeInfodf = RidgeInfo.to_pandas()
        
        LofarInfodf.to_csv(RLF.NLLR %source_name, columns = Lofarcolumns, header = True, index = False)
        RidgeInfodf.to_csv(RLF.NRLR %source_name, columns = Ridgecolumns, header = True, index = False)
        #np.savetxt('Ratios/NearestLofarLikelihoodRatios-%s.txt' %source_name, LofarInfo['AllWISE','Lofar_R_Distance', 'Lofar_LR', 'PossRA', 'PossDEC'], delimiter = ',', fmt='%s', encoding = 'utf-8')    
        #np.savetxt('Ratios/NearestRidgeLikelihoodRatios-%s.txt' %source_name, RidgeInfo['AllWISE','Ridge_R_Distance', 'Ridge_LR', 'PossRA', 'PossDEC'], delimiter = ',', fmt='%s', encoding = 'utf-8')   

#############################################

def LOFARLikelihoodRatio(table):
    """
    ***REDUNDANT***
    Calculates Best et. al.'s likelihood ratio using the R distances from 
    the LOFAR Catalogue position.
    NOTE: Can be used on the appropriate tables to calculate the 
    likelihood ratio for just the hosts.
    
    Parameters
    ----------
    
    table - array, shape(3,)
            the table of R distances and associated uncertainties on RA 
            and DEC.
    
    Returns
    -------
    
    likelihoodratio - list,
                      a list of likelihood ratios for the R distances from
                      the LOFAR Catalogue position.
    
    LRatDistance - array, shape(2,)
                   an array of the likelihood ratios and the corresponding
                   R distance from the LOFAR Catalogue position.
    
    """
    
    likelihoodratio = []
    distance = []
    
    for row in table:
        
        r, sigra, sigdec = RInfo(row)
        Y = Lambda(sigra, sigdec)
        
        LR = (np.float128(1.0 )/ (np.float128(2.0) * Y)) * np.exp(((r ** np.float128(2.0)) / np.float128(2.0)) * (np.float128(2.0) * Y - np.float128(1.0)))
        
        likelihoodratio.append(LR)
        distance.append(r)
    
    LRatDistance = np.column_stack((distance, likelihoodratio))
        
    return likelihoodratio, LRatDistance

#############################################
    
def NClosestDistances(available_sources, lofar_table, n):
    
    
    """
    Returns the number of hosts found within the Nth closest sources to the
    LOFAR catalogue position.  This R distance is  taken from the RDistances 
    file.  This function creates a .txt file containing the AllWise ID, 
    RDistance from Lofar Catalogue position, the calculated sigma's for RA 
    and DEC, the RA and DEC position for each possible source and the optical
    errs on RA and DEC for each of the n sources.
    
    Parameters
    ----------
    
    source_list - list,
                  a list of the sources that succesfully drew a ridgeline.         

    available_sources - 2D array, shape(x, 14), dtype='str'
                        array containing x source names and LOFAR catalogue
                        locations in pixel values, number of components
                        associated with source and LOFAR catalogue locations
                        as RA and DEC, the length of the source in pixels
                        and the total flux for all files available in a 
                        given directory. The host location in pixels, the
                        error on the RA and DEC of the host location, and
                        the RA and DEC of the host.
                        
    lofar_table - table,
                  the table containing all the LOFAR data from the .fits
                  catalogue.
                  
    n - int,
        the number of possible sources away from the LOFAR catalogue position,
        or the Nth closest source to the LOFAR catalogue position.

    Returns
    -------
    
    host_counter - float,
                   the number of hosts which are present within the Nth 
                   closest sources to the LOFAR catalogue position.
                  
    """
    
    #host_counter = 0
    #LofarRdistundern = []
    
    for asource in available_sources:
        source=asource['Source_Name']
        #ids = []
        LRdist = []
        sra = []
        sdec = []
        optra = []
        optdec = []
        raerr = []
        decerr = []
        #allwise = SourceInfo(source, lofar_table)[3]
        
        lmsize = asource['Size']/(3600.0*RLC.ddel) # pixels
                
        file5 = open(RLF.LDists %source, 'r')
        
        for row in file5:
            row_info = row.strip()
            row_cols = row_info.split(',')
            #allw = row_cols[0].strip("''").strip('b').strip("''")
            sRA = float(row_cols[0])
            sDEC = float(row_cols[1])
            distLR = float(row_cols[2])
            distLR /= lmsize
            optRA = float(row_cols[3])
            optDEC = float(row_cols[4])
            optRAerr = float(row_cols[5])
            optDECerr = float(row_cols[6])
            
            LRdist.append(distLR)
            sra.append(sRA)
            sdec.append(sDEC)
            optra.append(optRA)
            optdec.append(optDEC)
            raerr.append(optRAerr)
            decerr.append(optDECerr)
            #ids.append(allw)
            
        #sortids = [ids for LRdist, ids in sorted(zip(LRdist, ids))]
        sortra = [sra for LRdist, sra in sorted(zip(LRdist, sra))]
        sortdec = [sdec for LRdist, sdec in sorted(zip(LRdist, sdec))]
        sortoptra = [optra for LRdist, optra in sorted(zip(LRdist, optra))]
        sortoptdec = [optdec for LRdist, optdec in sorted(zip(LRdist, optdec))]
        sortraerr = [raerr for LRdist, raerr in sorted(zip(LRdist, raerr))]
        sortdecerr = [decerr for LRdist, decerr in sorted(zip(LRdist, decerr))]
        
        #nsortids = sortids[:n]
        nsortra = sortra[:n]
        nsortdec = sortdec[:n]
        nsortoptra = sortoptra[:n]
        nsortoptdec = sortoptdec[:n]
        nsortraerr = sortraerr[:n]
        nsortdecerr = sortdecerr[:n]
        nsortLRD = sorted(LRdist)[:n]
        rsiganderr = np.column_stack((nsortLRD, nsortra, nsortdec, nsortoptra, nsortoptdec, nsortraerr, nsortdecerr))
        
        columns = ['Lofar_R_Dist', 'Sigma_RA', 'Sigma_DEC', 'Opt_RA', 'Opt_DEC', 'Opt_RA_err', 'Opt_DEC_err']
        Nclosest = Table(rsiganderr, names = columns, dtype = ('f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))
        
        np.savetxt(RLF.NDist %source, Nclosest['Lofar_R_Dist','Sigma_RA', 'Sigma_DEC', 'Opt_RA', 'Opt_DEC', 'Opt_RA_err', 'Opt_DEC_err'], delimiter = ',', fmt='%s', encoding = 'utf-8')
        
        #for j in range(len(nsortids)):
        #    if allwise == nsortids[j]:
        #        host_counter += 1
                
        #for i in range(len(nsortLRD)):
        #    LofarRdistundern.append(nsortLRD[i])
    
        file5.seek(0)
    
    #return host_counter

#############################################
    
def NClosestSourcesDistances(source_list, lofar_table, n):
    
    """
    ***
    ISN'T NEEDED FOR DR2 AND WILL NOT WORK
    ***
    Returns the number of hosts found within the Nth closest sources to the
    LOFAR catalogue position. Also returns lists of the distances of the
    N closest sources from the LOFAR catalogue position, perpendicular 
    distance from the corresoponding ridgeline and distance along the 
    ridgeline.  This function uses the .txt file created in the the 
    CreateDistTable function.  Once these .txt files have been created 
    there is no need to run the previous function again.
    
    Parameters
    ----------
    
    source_list - list,
                  a list of the sources that succesfully drew a ridgeline.     
    
    lofar_table - table,
                  the table containing all the LOFAR data from the .fits
                  catalogue.    
    
    n - int,
        the number of possible sources away from the LOFAR catalogue position,
        or the Nth closest source to the LOFAR catalogue position.   

    Returns
    -------
    
    host_counter - float,
                   the number of hosts which are present within the Nth 
                   closest sources to the LOFAR catalogue position.
                 
    Ldistundern - list,
                  a list of all the distance from the LOFAR catalogue
                  position of all the sources which are the Nth closest to
                  to the LOFAR catalogue position.
                  
    Pdistundern - list,
                  a list of all the perpendicular distances from the ridgeline
                  of all possible sources which are the Nth closest to the 
                  LOFAR catalogue position.

    Adistundern - list,
                  a list of all the source distances along the ridgeline
                  of all possible sources which are the Nth closest to the 
                  LOFAR catalogue position.
                  
    """
    
    host_counter = 0
    Ldistundern = []
    Pdistundern = []
    Adistundern = []
    
    for source in source_list:
        ids = []
        Ldist = []
        Pdist = []
        Adist = []
        allwise = SourceInfo(source, lofar_table)[3]
        
        file4 = open(RLF.Dists %source, 'r')
        
        for row in file4:
            row_info = row.strip()
            row_cols = row_info.split(',')
            allw = row_cols[0].strip("''").strip('b').strip("''")
            distL = float(row_cols[2])
            distP = float(row_cols[3])
            distA = float(row_cols[4])
            
            Ldist.append(distL)
            Pdist.append(distP)
            Adist.append(distA)
            ids.append(allw)
            
        sortids = [ids for Ldist, ids in sorted(zip(Ldist, ids))]
        sortPD = [Pdist for Ldist, Pdist in sorted(zip(Ldist, Pdist))]
        sortAD = [Adist for Ldist, Adist in sorted(zip (Ldist, Adist))]
        
        nsortids = sortids[:n]
        nsortPD = sortPD[:n]
        nsortLD = sorted(Ldist)[:n]
        nsortAD = sortAD[:n]
        
        for j in range(len(nsortids)):
            if allwise == nsortids[j]:
                host_counter += 1
                
        for i in range(len(nsortLD)):
            Ldistundern.append(nsortLD[i])
            Pdistundern.append(nsortPD[i])
            Adistundern.append(nsortAD[i])
    
    file4.seek(0)
    
    return host_counter, Ldistundern, Pdistundern, Adistundern

#############################################
    
def PerpDistance(closest1, closest2, host):
    
    """
    Calculates the perpendicular distance from a point to a line.  Works
    out the perpendicular distance for a single host ID and it's ridgeline.
    If this distance lies outside the ridgeline then the shortest distance to
    the ridgeline is calculated from the closest point to the host.
    
    Parameters
    ----------
    
    host - array, shape(2, )
           the x and y pixel point of the host.
    
    m - float,
        the value of the gradient of the equation of the line joining
        the two closest points.
    
    c - float,
        the value of the intercept of the equation of the line joining
        the two closest points.
    
    Returns
    -------
    
    perpdis - float,
              the perpendicular distance of the host ID from a line.
    
    """
    
    Ix, Iy = PointOfIntersection(closest1, closest2, host)

    Ipix = np.array([Ix, Iy])
    perpdis = norm(Ipix - host)

    return perpdis
          
#############################################

def PointOfIntersection(closest1, closest2, host):
    
    """
    This function finds the point of intersection between the equation of the
    line joining the two points on the ridgeline closest to the host point
    and the equation of the perpendicular line through the host point.  If 
    this point lies outside the range of the two closest points, then the 
    distance from the two closest points is calculated and the nearest point
    is taken as the point of intercept.
    
    Parameters
    ----------
    
    closest1 - array, shape (2, ),
               the closest point to the host ID.
    
    closest2 - array, shape (2, ),
               the second closest point to the host ID.

    host - array, shape(2, )
           the x and y pixel point of the host.
    
    Returns
    -------
    
    Ix - float,
         the x pixel coordinate of the intercept.
         
    Iy - float,
         the y pixel coordinate of the intercept.
    
    """
##  To find a point of intersection on the RL to give the shorest distance (perpendicular) solve the
##  parametric vector equations.
##  Set up the rearranged answer and cooefficent matrices of the prametric equations for linalg.solve
##  from the rearrangement of the direction vectors.  It is perpendicular so the second direction is
##  given by (y, -x) from the first.
    Ans = np.array([host[0] - closest2[0], host[1] - closest2[1]])
    Cooefs = np.array([[closest1[0] - closest2[0], -(closest1[1] - closest2[1])],\
                           [closest1[1] - closest2[1], closest1[0] - closest2[0]]])
##  Solve the parametric equations for the values of t, s to give the point of intersection and 
##  calculate the point (only one needs to be done, so can remove one when I move in to the function).
    t, s = solve(Cooefs, Ans)
    Ix = host[0] + (closest1[1] - closest2[1]) * s
    Iy = host[1] - (closest1[0] - closest2[0]) * s
        
    if ((closest1[0] < Ix and Ix < closest2[0]) or (closest2[0] < Ix and Ix < closest1[0]))\
        and ((closest1[1] < Iy and Iy < closest2[1]) or (closest2[1] < Iy and Iy < closest1[1])):
        
        Ix = Ix  ##  This or the inequality needs reworking so it is a better piece of code!!
        Iy = Iy
        
    else:
        dis1 = norm(closest1 - host)
        dis2 = norm(closest2 - host)
        dis = min(dis1, dis2)
        
        if dis == dis1:
            Ix = closest1[0]
            Iy = closest1[1]
        else:
            Ix = closest2[0]
            Iy = closest2[1]
            
    return Ix, Iy

#############################################
        
def ProbHost(table):
    
    """
    ***REDUNDANT***
    Calculates the probability of finding a host between a distance of r and 
    r + dr as determined by Best et. al. (2003).
    Note:  It does not contain the dr factor as it is used later for the 
    division and the dr cancels out in this process.
    
    Parameter
    ---------
    
    table - array, shape(3,)
            the table of R distances and associated uncertainties on RA 
            and DEC.    
    
    Returns
    -------
    
    probhost - list,
               a list containing all the probabilties of there being a host
               lying between a distance r and r + dr.
    
    distprobhost - array, shape(2,)
                   an array containing all the probabilities of a host lying
                   between a distance r and r + dr away and the corresponding
                   r distance.    
    
    """
    
    probhost = []
    disthost = []
    
    for row in table:
        
        r, sigra, sigdec = RInfo(row)
        Y = Lambda(sigra, sigdec)
        
        prob = r * np.exp((-Y * (r ** 2)) / 2)
        
        probhost.append(prob)
        disthost.append(r)
    
    distprobhost = np.column_stack((disthost, probhost))
    
    return probhost, distprobhost

#############################################
    
def ProbPossSources(table):
    
    """
    ***REDUNDANT***
    Calculates the probabilty of a source lying between a distance of r and 
    r + dr as calculated by Best et. al (2003).
    Note:  It is missing the dr factor as it is later used for dividing and 
    the dr factor cancel out.
    
    Parameter
    ---------
    
    table - array, shape(3,)
            the table of R distances and associated uncertainties on RA 
            and DEC.
    
    Returns
    -------
    
    probposssources - list,
                      a list containing all the probabilties of there being a 
                      source lying between a distance r and r + dr.
    
    distprobposs - array, shape(2,)
                   an array containing all the probabilities of a source lying
                   between a distance r and r + dr away and the corresponding
                   r distance.
    
    """
    
    probposssources = []
    distposssources = []
    
    for row in table:
        
        r, SigRA, SigDEC = RInfo(row)
        Y = Lambda(SigRA, SigDEC)
        
        prob = 2 * Y * r * np.exp(-Y * (r ** 2))
        
        probposssources.append(prob)
        distposssources.append(r)
    
    distprobposs = np.column_stack((distposssources, probposssources))
        
    return probposssources, distprobposs

#############################################

def LDistance(DelRA, DelDEC):
    
    """
    Calculates the linear distance between two points using their given
    RA and DEC offsets.
    
    Parameters
    ----------
    
    ra - float,
         the RA in degrees of the position of interest on the sky (not used)

    dec - float,
          the Dec in degrees of the position of interest on the sky

    DelRA - float,
            the value of the offset of the RA coordinates between two points
            in arcseconds.
    
    DelDEC - float,
             the value of the offset of the DEC coordinates between two points
             in arcseconds.

    Returns
    -------

    l - float,
        the value of the distance between two given points, calculated using
        the RA and DEC offsets.
         
    """
    
    l = np.sqrt((DelRA ** 2) + (DelDEC ** 2))
    
    return l

#############################################

def RDistance(DelRA, DelDEC, SigRA, SigDEC):
    
    """
    Calculates the R distance from Best et. al. (2003) between two points
    using their given RA and DEC offsets and uncertainties.
    
    Parameters
    ----------
    
    DelRA - float,
            the value of the offset of the RA coordinates between two points
            in arcseconds.
    
    DelDEC - float,
             the value of the offset of the DEC coordinates between two points
             in arcseconds.    
    
    SigRA - float,
            the uncertainty on the RA in arcseconds.
            
    SigDEC - float,
             the uncertainty on the DEC in arcseconds.
    
    Returns
    -------

    r - float,
        the value of the distance between two given points, calculated using
        the RA and DEC offsets and uncertainties.
         
    """
    
    r = np.sqrt((DelRA ** 2 / SigRA ** 2) + (DelDEC ** 2 / SigDEC ** 2))
    
    return r

#############################################

def RInfo(row):
    
    """
    ***REDUNDANT***
    
    Retreives the information out of the tables created in HostRFromLofar
    and HostRFromRL for use in the likelihood ratio functions.
    
    Parameters
    ----------
    
    row - array, shape(3,),
          the row from the table containing the R distance and the 
          uncertainities.
    
    
    Returns
    -------
    
    r - float128,
        the value of the distance between two given points, calculated using
        the RA and DEC offsets and uncertainties.
    
    SigRA - float,
            the uncertainty on the RA in arcseconds.
            
    SigDEC - float,
             the uncertainty on the DEC in arcseconds.
             
    """
    
    r = row[0]
    SigRA = row[1]
    SigDEC = row[2]
        
    return np.float128(r), SigRA, SigDEC

#############################################
    
def SigmaR(LOFAR_RA_errRad, LOFAR_DEC_errRad, Host_RA_errOpt, Host_DEC_errOpt):
    
    """
    Calculates the uncertainities on the RA and DEC.  It uses the radio and
    optical uncertainties on the source and position away from it.  There is
    a preset assumption that the optical host has no radio erros and the radio
    point has no optical errors.  The uncertainties are in the same units as
    the inputted data.
    
    Parameters
    ----------
    
    LOFAR_RA_errRad - float,
                      the value of the radio uncertainity on the LOFAR
                      catalogue RA position in arcseconds.
            
    LOFAR_DEC_errRad - float,
                       the value of the radio uncertainity on the LOFAR
                       catalogue DEC posiiton in arcseconds.
                       
    Host_RA_errOpt - float,
                     the value of the optical uncertainity on the host RA
                     position in arcseonds.
                     
    Host_DEC_errOpt - float,
                      the value of the optical unceetainty on the host DEC
                      position in arcseconds.
    
    Constants
    ---------
    
    SigAst - float,
             the astrometric error from changing between the radio and 
             optical catalogues in arcseconds.
                  
    Returns
    -------
    
    SigRA - float,
            the uncertainty on the RA in arcseconds.
            
    SigDEC - float,
             the uncertainty on the DEC in arcseconds.
    
    """
    
    Host_RA_errRad = 0  ##  These are all currently set to 0 untilI have these values
    Host_DEC_errRad = 0  ##  or I know that these values are not needed because 
    LOFAR_RA_errOpt = 0  ##  they do not exist.
    LOFAR_DEC_errOpt = 0
    
    ##  NOTE: I know the square root and the square cancel each other however this is the
    ##  correct formula for working out the uncertainty and when creating the function I
    ##  kept forgetting that I had cancelled them out and thought I had worked the formula
    ##  out incorrectly.  In the hopes that I would stop doing that I have used the 
    ##  full formula.
    
    SigRA = np.sqrt((np.sqrt(LOFAR_RA_errRad ** 2 + Host_RA_errRad ** 2)) ** 2 + (np.sqrt(LOFAR_RA_errOpt ** 2 \
                + Host_RA_errOpt ** 2)) ** 2 + RLC.SigAst ** 2)
    
    SigDEC = np.sqrt((np.sqrt(LOFAR_DEC_errRad ** 2 + Host_DEC_errRad ** 2)) ** 2 + (np.sqrt(LOFAR_DEC_errOpt ** 2 \
                + Host_DEC_errOpt ** 2)) ** 2 + RLC.SigAst ** 2)

    return SigRA, SigDEC

#############################################
    
def SourceInfo(source, lofar_table):
    
    """
    Retrieves the LOFAR catalogue RA and Dec postion, AllWISE ID, ObjID
    redshift and size (in arcseconds) for the given source from the 
    lofar catalgoue table.
    
    Parameters
    ----------
    
    source - str,
             the LOFAR ID for the successfully drawn source.
                  
    lofar_table - table,
                  the table containing all the LOFAR data from the .fits
                  catalogue.
                  
    Returns
    -------
    
    lofara - float,
             the LOFAR catalogue RA postion.
             
    lofardec - float,
               the LFOA catalogue DEC postion.
    
    size - float,
           the size of the source in arcseconds.
    
    allwise - str,
              the AllWISE ID for the source.
    
    objid - float,
            the object ID number for the source.
    
    redz - float,
           the redshift for the source.
           
    """
    
    for row in lofar_table:
        if source == row[0]:
            lofarra = float(row[1])
            lofardec = float(row[2])
            #allwise = row[3]
            #objid = float(row[4])
            #redz = float(row[3])
            #size = float(row[4])   
    
    return lofarra, lofardec #size, redz

#############################################
    
def TableFromLofar(lofarcat):
    
    """
    Takes the required columns from the LOFAR catalogue of sources
    (all_eyeclass_4rms_v1.2_morph_agn_length_err.fits) and turns it into
    a table of data for use.  The columns in the catalogue are Source Name
    (the LOFAR Source Name), RA (the LOFAR RA position), DEC (the LOFAR DEC 
    position), AllWISE (the AllWISE ID), objID (the source object ID), z_best
    (the source redshift) and LGZ_size (the source size in arcseconds).

    Parameters
    -----------

    lofarcat - str,
               path and name to the catalogue .fits file.

    Returns
    -------
    
    Ltable - Table,
             the table containing all the LOFAR data from the .fits catalogue.
            
    """
    
    hdulist = fits.open(lofarcat)  ## Open the only Catalogue
    tbdata = hdulist[1].data  ## Find the data
    Name = tbdata.field(str(RLF.SSN))  ## Find the source name
    Lra = tbdata.field(str(RLF.SRA))  ## Find the LOFAR RA position column of the source (degree)
    Ldec = tbdata.field(str(RLF.SDEC))  ## Find the LOFAR DEC position column of the source (degree)
    #Lwise = tbdata.field('AllWISE')  ## Find the AllWISE ID
    #LID = tbdata.field('objID')  ## Find the object ID of the source
    #Lz = tbdata.field(str(RLF.LredZ))  ## Find the redshift of the source
    #Lsize = tbdata.field('LGZ_size')  ## Find the size of the source (arcsec)    
    source_names = np.column_stack((Name, Lra, Ldec))
    ## Stack all four columns next to each other. Note: Does it deal with missing data?
    columns = [str(RLF.SSN), str(RLF.SRA), str(RLF.SDEC)]  ## Creates column headings for calling rather than indices
    Ltable = Table(source_names, names = columns, dtype = ('S100', 'f8', 'f8'))  ## Turns it in to a table
    
    return Ltable

#############################################

def TableOfSources(catalogue):

    """
    Takes the catalogue of all PanStarrs and AllWISE sources
    (hetdex_ps1_allwise_radec.fits) and turns it into a table of data for
    use.  The columns in the catalogue are AllWISE (the AllWISE ID), objID
    (the source object ID), ra (the source RA position) and dec (the source
    DEC position).

    Parameters
    -----------

    catalogue - str,
                path and name to the catalogue .fits file.

    Returns
    -------
    
    table - Table,
            the table containing all the data from the .fits catalogue.
    

    """
    
    hdulist = fits.open(catalogue)  ## Open the only Catalogue
    tbdata = hdulist[1].data  ## Find the data
    WISE = tbdata.field(str(RLF.IDW))  ## Find the AllWISE ID
    UIDL = tbdata.field(str(RLF.IDP))  ## Find the object ID of the source
    #ID = tbdata.field(str(RLF.ID3))  ## Find the object ID of the source
    ra = tbdata.field(str(RLF.PossRA))  ## Find the RA position column of the source (degree)
    dec = tbdata.field(str(RLF.PossDEC))  ## Find the DEC position column of the source (degree)
    #ra_err = RLF.PRAErr  ## Find the error on the RA position (arcsec)
    #dec_err = RLF.PDECErr  ## Find the error on the DEC position (arcsec)
    w1 = tbdata.field(str(RLF.OptMagA))  ## Find the W1 magnitude
    #w1_err = tbdata.field('W1magErr')  ## Find the error on the W1 magnitude
    r = tbdata.field(str(RLF.OptMagP))  ## Find the i magnitude
    #i_err = tbdata.field('iErr')  ## Find the error on the i magnitude
    source_names = np.column_stack((WISE, UIDL, ra, dec, w1, r))
    ## Stack all four columns next to each other. Note: Does it deal with missing data?
    columns = [str(RLF.IDW), str(RLF.IDP), str(RLF.PossRA), str(RLF.PossDEC), str(RLF.OptMagA), str(RLF.OptMagP)]  ## Creates column headings for calling rather than indices
    table = Table(source_names, names = columns, dtype = ('S100', 'S100', 'f8', 'f8', 'f8', 'f8'))  ## Turns it in to a table
    
    return table

#############################################

"""
Reference
---------
Best et. al. (2003), CENSORS: A Combined EIS-NVSS Survey Of Radio Sources.
I. Sample definition, radio data and optical identifications, MNRAS,
346, 627 - 683

"""

#############################################

