# ridge_toolkit.py
# Toolkit for ridge-line finding in LOF#R extended sources
# Based on Jude's original code, modified by Beatriz
# Further modified by Joanna
# Modified some more by Bonny

import RidgelineFilesDR2 as RLF
import RLConstantsDR2 as RLC
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.wcs import WCS, utils
from math import atan2
from numpy import pi, cos, sin, sqrt, absolute
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pyregion
from scipy import ndimage
from skimage import img_as_uint
from skimage.measure import label
from skimage.morphology import octagon
from skimage.morphology import erosion
from skimage.feature import peak_local_max
from skimage.filters import threshold_minimum
import copy
from sizeflux_tools import Flood,length3,Mask

#############################################

def AngleRange(phi, dphi_plus, dphi_minus):

    """
    Calculates the upper and lower angle limits in [-pi,pi]
    range, given the mid-angle value and deviation allowance
    in both directions

    Parameters
    -----------

    phi - float,
          angle, from which the range is calculated, given in radians

    dphi_plus - float,
                deviation from phi in the direction of increasing phi,
                given in radians

    dphi_minus - float,
                 deviation from phi in the direction of decreasing phi,
                 given in radians

    Returns
    -----------

    phirange - 1D array, shape(2, )
               array of angle limits in the [-pi,pi] range in the form
               [lower_value, upper_value]
    """


    phi_up = PiRange(phi + dphi_plus)
    phi_down = PiRange(phi - dphi_minus)

    phirange = np.array([phi_down, phi_up])

    return phirange

#############################################

def AreaFluxes(fluxes):

    """
    Convolves the original 2D array of flux values at single
    pixels with a centre-heavy cross kernel.

    Parameters
    -----------

    fluxes - 2D array,
             array containing flux values for single pixels

    Constants
    ---------
    
    KerW - int,
           The weight of the kernel in the convolution

    Returns
    -----------

    area_fluxes - 2D array,
                  input array convolved with a centre-heavy cross
                  kernel

    Notes
    -----------

    Uses 'constant' image boundaries with the values outside of the
    original array set to 0.0

    """
    fluxes = np.ma.filled(fluxes, 0)
    fluxes_data = np.nan_to_num(fluxes)
    kernel = ndimage.generate_binary_structure(2,1).astype(np.int)
    kernel[1,1] = RLC.KerW    
    area_flux = ndimage.convolve(fluxes_data, kernel, mode='constant', cval=0.0)
    area_flux = area_flux/(RLC.KerW + 4)
    area_flux[area_flux == 0] = np.nan
    
    return area_flux

#############################################

def CheckQuadrant(angle):

    """
    Returns angle quadrant for given input angle

    Parameters
    -----------

    angle - float,
            angle of interest, given in radians

    Returns
    -----------

    quadrant - int,
               quadrant, in which the angle is located.
    """

    if 0 <= angle <= pi/2.:
        quadrant = 1

    elif pi/2. <= angle <= pi:
        quadrant = 2

    elif -pi <= angle <= -pi/2.:
        quadrant = 3

    elif -pi/2. <= angle <= 0:
        quadrant = 4

    return quadrant

#############################################

def ComponentsTable(components):
    
    """Returns a table of LOFAR Galaxy Zoo Associated components to each 
    source in available sources, their central position RA and DEC, 
    their Total Flux, the Maj and Min axes and the angle of rotation of
    the ellipse.

    This is converted from the standard component table by putting the
    fields in a specific order so that they can be indexed numerically.
    
    Parameters
    ----------
    
    components - str,
                 path and name to the components catalogue .fits file
    
    Returns
    -------
    
    CompTable - astropy table, shape ( , 8),
                A table containing all the components linked to each 
                source.  Listing the location RA and DEC, the ellipse
                Maj and Min axes and the angle of rotation.  As well as
                the Total Flux of the component.

    """
    
    t = Table.read(components)
    columns2 = [RLF.CSN, RLF.CCN, RLF.CRA, RLF.CDEC, RLF.CTF, RLF.CMAJ, RLF.CMIN, RLF.CPA]
    component_list = [t[c] for c in columns2]
    CompTable = Table(component_list, names = columns2)
    
    return CompTable
    
#, dtype = ('s100', 's100', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8') goes after names = columns2 incase I want to put it back in

#############################################

def ConfigurePlots():

    """
    
    Setting global Matplotlib settings
    
    Constants
    ---------
    
    lw - int,
         linewidth of the graph
    col - str,
          line color
    ms - int,
         markersize
    mt - str,
         mathtext font
    ff - str,
         the choice of font family
    fs - int,
         the font size
    ts - int,
         the size of the title
    ls - int,
         the size of the labels
    ym - float,
         the sixe of the y-margin
    figs - float,
           the size of the figure, will be a square if inserted twice
    fonts - str,
            the font size of the legend, e.g. small, medium, large
    sav - str,
          save format (default is eps)
    
    
    """
    
    rcParams['lines.linewidth'] = RLC.lw
    rcParams['lines.color'] = RLC.col
    rcParams['lines.markersize'] = RLC.ms
    rcParams['mathtext.rm'] = RLC.mt
    rcParams['font.family'] = RLC.ff
    rcParams['font.size'] = RLC.fs
    rcParams['axes.titlesize'] = RLC.ts
    rcParams['axes.labelsize'] = RLC.ls
    rcParams['axes.ymargin'] = RLC.ym
    rcParams['figure.figsize'] = RLC.figs, RLC.figs
    rcParams['legend.fontsize'] = RLC.fonts
    rcParams['savefig.format'] = RLC.sav
    rcParams['pcolor.shading'] = RLC.pcol
    
#############################################   

def CreateCutouts(available_sources):
    
    """
    This function creates smaller sized cutouts centred on the LOFAR
    catalogue position and out to lmsize from Mingo et al. in all
    directions for the FITS and RMS arrays.  It allows the program to
    run on a smaller cutout making it run quicker and eliminating a
    majority of excess sources before the flood fill and masking steps
    take place, allowing them to run quicker.  This function only needs
    to be run once, as once the cutouts are generated this does not
    need to be duplicated.
    
    Parameters
    ----------
    
    available_sources - astropy Table
    
    Constants
    ---------
    
    Start - int,
            the index to start the selection of sources from (default = 0)
    
    Step - int,
           the fraction step size of the function. 1 will run everything,
           2 half etc
    
    ImSize - float,
             the fraction of the source size that the final image is cut
             down to. 1 will match the size of the source from the centre
             of the cutout
    
    Returns
    -------
    
    Saves the smaller cutouts to the fits_cutouts and rms4_cutouts folders
    which are then accessed through the program.
    
    """
    
    for source in available_sources[RLC.Start::RLC.Step]:
        
        source_name = source['Source_Name']
        centre_pos = SkyCoord(source['RA']*u.deg,source['DEC']*u.deg,frame='fk5')
        lmsize = source['Size']/(3600.0*RLC.ddel) # convert to pixels
        print('Making cutout for source',source_name,'with size',lmsize,'pixels')
        print('    centre_pos is',centre_pos)
        hdu = DefineHDU(source_name)
        data = hdu[0].data
        size = (2 * lmsize, 2 * lmsize)
        wcs = WCS(hdu[0].header)  ## Keep world coordinate system
        cutout = Cutout2D(data, centre_pos, size, wcs = wcs)
        hdu[0].data = cutout.data
        hdu[0].header.update(cutout.wcs.to_header())
        hdu[0].writeto(str(RLF.fitsfile) + source_name + '-cutout.fits', overwrite = True)
        
        flux_array = GetFluxArray(source_name)
        cutout2 = Cutout2D(flux_array, centre_pos, size, wcs=wcs)
        np.save(str(RLF.npyfile) + source_name + '-cutout.npy', cutout2.data)
        

#############################################

def DefineCutoutHDU(source_name):

    """
    Returns an opened .fits file cutout from a given source loaded 
    from an available .fits file

    Parameters
    -----------
    source_name - str,
                  source name used to identify a .fits file


    Returns
    -----------
    hduC_open - an opened .fits file

    """

    cutout_name = str(RLF.fitscutout) + source_name + '-cutout.fits'
    hduC_open = fits.open(cutout_name)


    return hduC_open    

############################################# 

def DefineHDU(source_name):

    """
    Returns an opened .fits file from a given source loaded from an
    available .fits file

    Parameters
    -----------
    source_name - str,
                  source name used to identify a .fits file


    Returns
    -----------
    hdu_open - an opened .fits file

    """

    fits_name = str(RLF.fits) + source_name + '.fits'
    hdu_open = fits.open(fits_name)

    return hdu_open

#############################################

def ErodedMaxima(area_fluxes):

    """
    Returns a numpy array with local maxima locations found from eroding
    the original image as well as the the eroded structure, to be used
    in weighted angle averaging

    Parameters
    -----------

    area_fluxes - 2D array,
                  raw image array convolved with a centre-heavy cross
                  kernel

    Constants
    ---------
    
    TRel - float,
           the value that relative threshold must be greater than in 
           the erosion
    pThresh - float,
              peak threahold value for the erosion
    mindis - float,
             minimum distance between peaks in erosion (5 = beam size)
    Octm - int,
           the m parameter of octagon used to erode, defines the size 
           of the horizontal and vertical sides
    Octn - int,
           the n parameter of octagon used to erode, defines the size
           of the height and width of the slanted sides
    
    Returns
    -------

    maxima - 2D array, shape(n, 2)
             array of n peak positions in form [row, col],
             which corresponds to [y-pix, x-pix] when presented on a graph

    maxima_array - 2D array, shape(area_fluxes.shape)
                   array of boolean values with pixels containing
                   a maximum having value = True

    eroded - 2D array,
             input array transformed into dtype='uint16', eroded using
             an octagonal structuring element of size 9x9 pixels

    """


    zero_filled = np.ma.filled(np.ma.masked_invalid(area_fluxes), 0)
    
    try:
        thresh = threshold_minimum(zero_filled)
    except RuntimeError:
        print('Runtime Error')
    else:
        thresh_rel = thresh/np.nanmax(area_fluxes)
        peak_thresh = RLC.pThresh
        min_dist = RLC.mindis

        area_intfluxes = img_as_uint(area_fluxes)
        eroded = erosion(area_intfluxes, octagon(RLC.Octm,RLC.Octn))

    # include localised high intensity peaks, which dominate
    # low-intensity background and get eroded away

        if thresh_rel > RLC.TRel:

            cutoff = np.ma.filled(np.ma.masked_where(zero_filled < thresh,\
                    area_intfluxes), 0)
            eroded = np.ma.filled(np.ma.masked_where(zero_filled > thresh,\
                    eroded), 0)
            eroded += cutoff
            peak_thresh *= thresh_rel

        maxima_original = peak_local_max(zero_filled, min_distance=min_dist,\
                threshold_rel=peak_thresh)
        maxima_eroded = peak_local_max(eroded, min_distance=min_dist, \
                threshold_rel=peak_thresh)

        maxima = maxima_eroded
        maxima_array = np.full_like(area_fluxes, False, dtype=bool)

        if maxima.ndim == 1:
            if maxima_original.ndim == 1:
                diff = maxima_eroded - maxima_original
                diff_r = sqrt(np.power(diff[0], 2) + np.power(diff[1], 2))

                if diff_r <= 5.:
                    maxima = maxima_original
            else:
                diff_array = maxima_eroded - maxima_original
                diff_r = sqrt(np.power(diff_array[:,0], 2) + \
                        np.power(diff_array[:,1], 2))

                if np.amin(diff_r) <= 5.:
                    maxima = maxima_original[np.argmin(diff_r),:]

            maxima_array[int(maxima[0]), int(maxima[1])] = True

        else:
            for idx in range(maxima_eroded.shape[0]):

                diff_array = maxima_eroded[idx,:] - maxima_original
                diff_r = sqrt(np.power(diff_array[:,0], 2) + \
                        np.power(diff_array[:,1], 2))

                if np.amin(diff_r) <= 5.:
                    maxima[idx,:] = maxima_original[np.argmin(diff_r),:]

                maxima_array[int(maxima[idx,0]), int(maxima[idx,1])] = True

        maxima = maxima.astype('float64')
        maxima += 0.5
        maxima_unique = np.unique(maxima, axis=0)

    #eliminate accidental maxima islands, which escape the image processing algorithm

        idxmax = maxima_unique.shape[0]
        idx = 0
        while idx < idxmax:

            if idx == idxmax - 1:
                idx += 1

            else:
                dx = int(absolute(maxima_unique[idx,1]-maxima_unique[idx+1,1]))
                dy = int(absolute(maxima_unique[idx,0]-maxima_unique[idx+1,0]))

                if (dx == 1 and dy == 0)\
                    or (dx == 0 and dy == 1) \
                    or (dx == 1 and dy == 1):

                    maxima_array[int(maxima_unique[idx,0] - 0.5), \
                            int(maxima_unique[idx,1] - 0.5)] = False
                    maxima_unique = np.delete(maxima_unique, (idx), axis=0)
                    idxmax -= 1

                else:
                    idx += 1

        maxima_out = maxima_unique

        return maxima_out, maxima_array, eroded

#############################################

def FindRidges(area_fluxes, init_point, R, dphi, lmsize, \
               CompTable, n_comp, source_name, hdu):

    """
    Returns numpy arrays with point coordinates along 2 ridge lines
    and seperate numpy arrays containing the angular directions and
    length along the ridgeline associated with these points. The use
    of Rtot in the loops allows for the summation of the step sizes
    to calculate the length of the ridge lines in pixels.

    Parameters
    -----------

    area_fluxes - 2D array,
                  raw image array convolved with a centre-heavy cross
                  kernel
    
    init_point - 1D array, shape(2,)
                 starting point position for ridge detection    
              
    R - float,
        ridge finding step, given in pixels

    dphi - float,
           half of the value of the cone opening angle
    
    lmsize - str,
             the size of the source in pixels
    
    Constants
    ---------
    
    Lcone - float,
            the angle in degrees of the larger search cone
    MaxLen - float,
             the multiplier of the source length (lmsize) to determine
             the maximum ridgeline length
    ipit - int,
           the number of interations of Retry Directions and initial
           point finding that the code should attempt before determining
           a problematic source due to initial point finding

    Returns
    -------

    ridge1 - 2D array,
             array of ridge point coordinate values found
             for direction 1, in the form of [x_coord, y_coord]

    phival1 - 1D array,
              array of angular directions associated with each
              ridge point found for direction 1. Values in the
              range [-pi,pi]
    
    Rlen1 - 1D array,
            array of lengths associated with each ridge point
            found for direction 1.

    ridge2 - 2D array,
             array of ridge point coordinate values found
             for direction 2, in the form of [x_coord, y_coord]

    phival2 - 1D array,
              array of angular directions associated with each
              ridge point found for direction 2. Values in the
              range [-pi,pi]
              
    Rlen2 - 1D array,
            array of lengths associated with each ridge point
            found for direction 2.
            
    Notes
    -----------
    Returns np.nan arrays for both ridge1 and ridge2 if unable
    to find #ial directions for ridge search.

    The phival1 and phival2 arrays can be directly used to decide
    about 'bendiness' of the source.
    The ridge point locations can be easily turned into pixel labels
    within pixel array by using np.floor(ridge1).astype('int') command.
    This allows for direct extraction of an intensity profile along
    the ridgeline

    """
    
    #print(float(lmsize))
    #print(0.75 * float(lmsize))
    
    Error = 'N/A'
    try:
        maxima, maxima_array, eroded = ErodedMaxima(area_fluxes)
    except:
        print('Erosion Error')
        Error = 'Erosion'
    else:
    # leaves the option of working with the eroded image only.
    # It hasn't proved useful for the HETDEX data set
        new_fluxes = area_fluxes
        #new_fluxes = np.ma.masked_array(area_fluxes, \
        #mask=np.ma.masked_equal(eroded, 0).mask, copy=True)
        if RLC.debug == True:
            print('Initial Cones')
        phi_val1, phi_val2, cone1, cone2, init_point, Error = InitCones(area_fluxes,\
                                init_point, np.radians(75), lmsize, CompTable, n_comp, source_name, hdu)

        #phitot1 = dphi + phi_val1  ##part of the angle restriction idea
        #phitot2 = dphi + phi_val2  ##part of the angle restriction idea
        if RLC.debug == True:
            print('Initial Cones Completed')
        #print(Error, phi_val1, phi_val2)
        
        if (np.isnan(phi_val1) and np.isnan(phi_val2)):

            ridge1 = ridge2 = np.array([np.nan, np.nan])
            Rlen1 = np.array([np.nan, np.nan])
            Rlen2 = np.array([np.nan, np.nan])

        else:
            larger_cone1, larger_cone2, init_point, Error = InitCones(area_fluxes, \
                            init_point, np.radians(RLC.Lcone), lmsize, CompTable, n_comp, source_name, hdu)[-4:]
            
            if RLC.debug == True:
                print('Larger Init Cones completed')
            
            if type(larger_cone1) == np.ndarray or type(larger_cone2) == np.ndarray:

                ridge1 = ridge2 = np.array([np.nan, np.nan])
                Error = 'Unable_to_Find_Initial_Directions'
                Rlen1 = np.array([np.nan, np.nan])
                Rlen2 = np.array([np.nan, np.nan])
                
            else:
                
                chain_mask1 = larger_cone1.mask
                chain_mask2 = larger_cone2.mask
    
                new_point1, new_phi1, chain_mask1, RFNew1 = GetFirstPoint(new_fluxes, \
                                                                  init_point, cone1, chain_mask1, R)
                ridge1 = np.vstack((init_point, new_point1))
                new_point2, new_phi2, chain_mask2, RFNew2 = GetFirstPoint(new_fluxes, \
                                                                  init_point, cone2, chain_mask2, R)
                ridge2 = np.vstack((init_point, new_point2))
                
                #print('RFNew1 = ' + str(RFNew1))
                #print('RFNew2 = ' + str(RFNew2))
                
                #print('Larger Cones Completed')
    
                if np.isnan(new_phi1) or np.isnan(new_phi2):
    
                    print('Unable to find the first ridge point. '\
                          'Further source analysis is aborted.')
                    ridge1 = ridge2 = np.array([np.nan,np.nan])
                    #ridge2 = np.full_like(init_point, np.nan)
                    Error = 'Unable_to_Find_First_Ridge_Point'
                    Rlen1 = np.array([np.nan, np.nan])
                    Rlen2 = np.array([np.nan, np.nan])                
    
                else:
                    
                    phi_val1 = np.array([0, new_phi1])
                    phi_val2 = np.array([0, new_phi2])
                    Rmax = RLC.MaxLen * float(lmsize)
    
                    if RLC.debug == True:
                        print('Tracing first ridgeline')
                    #Rcounter1 = 0
                    Rtot1 = RFNew1
                    Rlen1 = np.array([0, Rtot1])
                
                    while True and Rtot1 < Rmax:# and Rcounter1 <= 15: #and phitot1 < 24 * dphi: ##Trying to add an angle restriction
                        #Rtot1 += R
                        #Rcounter1 += 1
                        cone1, check_cone1 = GetRidgeCone(area_fluxes, \
                                                          maxima_array, ridge1[-1,:], \
                                                          new_phi1, chain_mask1, dphi)
                        new_point1, new_phi1, chain_mask1, RNew1 = GetRidgePoint(new_fluxes, ridge1[-1,:], \
                                                                  cone1, check_cone1, chain_mask1, R, lmsize, Rtot1, Rmax)
                        Rnewlen1 = Rtot1 + RNew1
                        Rlen1 = np.append(Rlen1, Rnewlen1)
                        ridge1 = np.vstack((ridge1, new_point1))
                        phi_val1 = np.append(phi_val1, new_phi1)
                        #print('RNew1 = ' + str(RNew1))
                        Rtot1 += RNew1
    
                        #phitot1 += float(phi_val1) ## Part of trying to add an angle restriction
                        
                        if np.isnan(new_phi1):
                            break
                    
                    if RLC.debug == True:
                        print('Tracing second ridgeline')
                    Rtot2 = RFNew2
                    Rlen2 = np.array([0, Rtot2])
                    #Rcounter2 = 0
                    while True and Rtot2 < Rmax:# and Rcounter2 <= 15: #and phitot2 < 24 * dphi: ##Trying to add an angle restriction
                        #Rtot2 += R
                        #Rcounter2 += 1
                        cone2, check_cone2 = GetRidgeCone(area_fluxes, maxima_array, ridge2[-1,:], \
                                                          new_phi2, chain_mask2, dphi)
                        new_point2, new_phi2, chain_mask2, RNew2 = GetRidgePoint(new_fluxes, ridge2[-1,:], \
                                                                  cone2, check_cone2, chain_mask2, R, lmsize, Rtot2, Rmax)
                        Rnewlen2 = Rtot2 + RNew2
                        Rlen2 = np.append(Rlen2, Rnewlen2)
                        ridge2 = np.vstack((ridge2, new_point2))
                        phi_val2 = np.append(phi_val2, new_phi2)
                        #print('RNew2 = ' + str(RNew2))
                        Rtot2 += RNew2                    
                        
                        #phitot2 += float(phi_val2) ##Part of trying to add an angle restriction
                        
                        if np.isnan(new_phi2):
                            break
                
                #print('Rtot1 = ' + str(Rtot1))
                #print('Rtot2 = ' + str(Rtot2))
        return ridge1, phi_val1, Rlen1, ridge2, phi_val2, Rlen2, Error

#############################################


#############################################

def GetAvailableSources(filename):

    """
    Returns a numpy array with source names, LOFAR catalogue locations 
    in pixel values, number of components, LOFAR catalogue locations RA 
    and DEC, size of source in pixels, and Host ID location in pixels for
    available sources. It uses the fits file header and the region file 
    to convert RA, DEC ID positions into pixel locations and source size
    in pixel length.

    Parameters
    -----------
    filename - str,
               full name of the file containing a list of sources and
               optical ID locations in RA and DEC coordinates and the
               length of the source in arcsecs
    
    Constants
    ---------
    
    rdel - float,
           the equivalent pixel value to RA in a FITS file, set value
    
    ddel - float,
           the equivalent pixel value to DEC in a FITS file, set value

    Returns
    -----------
    available_sources - 2D array, shape(n, 14), dtype='str'
                        array containing n source names and LOFAR catalogue
                        locations in pixel values, number of components
                        associated with source and LOFAR catalogue locations
                        as RA and DEC, the length of the source in pixels
                        and the total flux for all files available in a 
                        given directory. The host location in pixels, the
                        error on the RA and DEC of the host location, and
                        the RA and DEC of the host.

    Notes
    -----------
    The output array is of dtype='str', hence the pixel locations of the
    optical IDs are cast as strings. Care must be taken in reading these
    values later in the processing.  I.E. Remember the lmsize needs float.

    """

    sources_list = open(filename,'r')
    available_sources = np.empty(10)

    for source in sources_list:

        source_info = source.strip()
        info_cols = source_info.split()
        source_name = info_cols[0].strip("''").strip('b').strip("''")
        source_name = source_name.rstrip()
        Lrx = float(info_cols[1])
        Lry = float(info_cols[2])
        ncomp = info_cols[3]
        size = float(info_cols[4])
        flux = float(info_cols[5])
        #idrx = float(info_cols[6])
        #idry = float(info_cols[7])
        LraE = float(info_cols[6])
        LdecE = float(info_cols[7])
        fits_name = str(RLF.fits) + source_name + '.fits'
        numpy_name = str(RLF.rms4) + source_name + '.npy'


        if os.path.isfile(fits_name) and os.path.isfile(numpy_name):

            hdu = fits.open(fits_name)
            header = hdu[0].header
            wcs = WCS(header)
#            ra = hdu[0].header['CRVAL1']
#            dec = hdu[0].header['CRVAL2']
#            xra = hdu[0].header['CRPIX1']
#            yra = hdu[0].header['CRPIX2']
#            data = hdu[0].data
#            ymax,xmax = data.shape
#            xr = np.array(range(xmax))
#            yr = np.array(range(ymax))
#            xg,yg = np.meshgrid(xr,yr)
#            xmin,xmax = xg.min(),xg.max()
#            ymin,ymax = yg.min(),yg.max()

            # convert RA, DEC coordinates into pixel locations for LOFAR
#            coords = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, \
#                    frame='fk5')
#            coords2 = SkyCoord(ra=float(Lrx)*u.degree, \
#                    dec=float(Lry)*u.degree, frame='fk5')
#            dist_x, dist_y = coords.spherical_offsets_to(coords2)
#            dx = float((dist_x.to(u.arcsec)/u.arcsec)/(3600.0*RLC.rdel))
#            dy = float((dist_y.to(u.arcsec)/u.arcsec)/(3600.0*RLC.ddel))
#            source_Lxpix = xra + int(dx) - (xmin-1)
#            source_Lypix = yra + int(dy) - (ymin-1)
#
#            source_Lxpix += 0.5
#            source_Lypix += 0.5
#            
#            # convert RA, DEC coordinates into pixel locations for Optical ID
#            coords3 = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, \
#                    frame='fk5')
#            coords4 = SkyCoord(ra=float(idrx)*u.degree, \
#                    dec=float(idry)*u.degree, frame='fk5')
#            dist_idx, dist_idy = coords3.spherical_offsets_to(coords4)
#            iddx = float((dist_idx.to(u.arcsec)/u.arcsec)/(3600.0*RLC.rdel))
#            iddy = float((dist_idy.to(u.arcsec)/u.arcsec)/(3600.0*RLC.ddel))
#            source_idxpix = xra + int(iddx) - (xmin-1)
#            source_idypix = yra + int(iddy) - (ymin-1)
#
#            source_idxpix += 0.5
#            source_idypix += 0.5
            
            #hostcoords = SkyCoord(float(idrx) * u.degree, float(idry) * u.degree, \
                                  #frame = 'fk5')
            #hostpix = utils.skycoord_to_pixel(hostcoords, wcs, origin = 1)
            #source_idxpix = hostpix[0]
            #source_idypix = hostpix[1]

            lofarcoords = SkyCoord(float(Lrx) * u.degree, float(Lry) * u.degree, \
                                  frame = 'fk5')
            lofarpix = utils.skycoord_to_pixel(lofarcoords, wcs, origin = 1)
            source_Lxpix = lofarpix[0]
            source_Lypix = lofarpix[1]
            
            lengthpix = (size/RLC.ddel/3600.0)  ## Converting the length of the source in to pixels

            source_summary = np.array([source_name, source_Lxpix, \
                    source_Lypix, ncomp, Lrx, Lry, lengthpix, flux, \
                    LraE, LdecE])
            available_sources = np.vstack((available_sources, \
                    source_summary))

    available_sources = available_sources[1:,:]

    return available_sources

#############################################

def GetCone(area_fluxes, coord_origin, phirange):

    """
    Returns a search cone based on the search point coordinates
    and the range of angle values.

    Parameters
    -----------
    area_fluxes - 2D array,
                  raw image array convolved with a centre-heavy cross
                  kernel

    coord_origin - 1D array, shape(2, )
                   cone origin point

    phirange - 1D array, shape(2, )
               array of angle limits in the [-pi,pi] range in the form
               [lower_value, upper_value]

    Returns
    -----------

    cone - 2D masked array,
           area_fluxes array masked to preserve the search cone only

    """

    r, phi = PolarCoordinates(area_fluxes, coord_origin)

    if 3 <= CheckQuadrant(phirange[1]) <= 4 and \
    1 <= CheckQuadrant(phirange[0]) <= 2:
        phi = np.ma.masked_inside(phi, phirange[1], phirange[0])

    else:
        phi = np.ma.masked_outside(phi, phirange[0], phirange[1])

    cone = np.ma.masked_array(area_fluxes, mask=phi.mask)

    return cone

#############################################

def GetCutoutArray(source_name):

    """
    Returns a cutout numpy array for a given source loaded 
    from an available .npy file

    Parameters
    -----------
    source_name - str,
                  source name used to identify a .npy file


    Returns
    -----------
    cutout_array - 2D array,
                   array of pixel fluxes loaded from file

    """

    Cnumpy_name = str(RLF.rms4cutout) + source_name + '-cutout.npy'
    cutout_array = np.load(Cnumpy_name)

    return cutout_array    

#############################################

def GetFirstPoint(area_fluxes, init_point, cone, chain_mask, R):

    """
    Returns a first point along the ridge with its angular direction
    and updated ridge chain mask (R pix radius mask around ridge
    points found in given direction). It returns the length of the
    ridge step for summation in FindRidges.

    Parameters
    -----------
    area_fluxes - 2D array,
                  raw image array convolved with a centre-heavy cross
                  kernel

    init_point - 1D array, shape(2, )
                 initial ridge finding point

    cone - 2D masked array,
           area_fluxes array masked to preserve the search cone only

    chain_mask - masked array mask,
                 already existing mask, which covers an R pix radius
                 around previous ridge points. In place to avoid the
                 ridge line wrapping onto itself

    R - float,
        ridge finding step, given in pixels


    Returns
    -----------

    new_point - 1D array, shape(2, )
                first point along the ridge line

    new_phi - float,
              angular direction of the new ridge point, found by the
              weighted average of polar coordinate angles at the given
              search point

    chain_mask - masked array mask,
                 updated input mask
                 
    RFNew - float,
            length of first ridge step, given in pixels

    """

    r, phi = PolarCoordinates(area_fluxes, init_point)
    Rmax = np.ma.max(r)
    cone_mask = cone.mask
    r_mask = np.ma.mask_or(np.ma.masked_outside(r, 0, R).mask, \
                           chain_mask)
    slice_mask = np.ma.mask_or(r_mask, cone_mask)
    cone_slice = np.ma.masked_array(area_fluxes, mask=slice_mask, \
                                    copy=True)

    safety_stop = False
    RFNew = 0

    while np.ma.count(cone_slice) == 0:

        R += 1
        r_mask = np.ma.mask_or(np.ma.masked_outside(r, 0, R).mask, \
                               chain_mask)
        slice_mask = np.ma.mask_or(r_mask, cone_mask)
        cone_slice = np.ma.masked_array(area_fluxes, mask=slice_mask, \
                                        copy=True)
        RFNew = R
        

        if R > Rmax:
            safety_stop = True
            new_phi = np.nan
            break

    if safety_stop == True:

        new_point = np.array([np.nan, np.nan])

    else:
        cone_phislice = np.ma.masked_array(phi, mask=slice_mask, \
                                           copy=True)
        quad1, quad2 = CheckQuadrant(np.ma.min(cone_phislice)), \
                        CheckQuadrant(np.ma.max(cone_phislice))
        diff = np.ma.max(cone_phislice) - np.ma.min(cone_phislice)

        if 3 <= quad1 <= 4 and 1 <= quad2 <= 2 and diff > pi:

            for idx, phi in np.ndenumerate(cone_phislice):
                if cone_phislice[idx] is np.ma.masked:
                    continue

                else:
                    phi = TwoPiRange(phi)
                    cone_phislice[idx] = phi

        cone_flux = np.ma.sum(cone_slice)
        weighted_phi = np.sum(np.multiply(cone_phislice, cone_slice))
        average_phi = PiRange(weighted_phi/cone_flux)
        new_point = np.array([R*cos(average_phi), R*sin(average_phi)]) + init_point
        new_phi = average_phi

        chain_mask = np.ma.mask_or(np.ma.masked_inside(r, 0, R).mask, chain_mask)
        RFNew = R
      

    return new_point, new_phi, chain_mask, RFNew

#############################################

def GetFluxArray(source_name):

    """
    Returns numpy array for a given source loaded from an
    available .npy file

    Parameters
    -----------
    source_name - str,
                  source name used to identify a .npy file


    Returns
    -----------
    flux_array - 2D array,
                 array of pixel fluxes loaded from file

    """

    numpy_name = str(RLF.rms4) + source_name + '.npy'
    flux_array = np.load(numpy_name)

    return flux_array

#############################################

def GetMaskedComp(hdu, source, components, flux_array, CompTable):
    
    """
    Masks out unrelated components.  Searches for components within a cutout 
    distance of the source and masks them. Any components outside of the distance
    are not included in the masking process.  
    It creates a table from the component catalogue of Component Name, RA, DEC,
    Total Flux, Major and Minor Axes and PA (the angle of rotation) of the ellipse 
    associated with each component. 
    Finds the components associated to the source and then searches and masks 
    the components not associated to the source.
    
    Parameters
    ----------
    
    source - row from available_sources, len = 6
             contains source name, optical poistion in x and y pixels, RA and 
             DEC and
             number of components associated to source.
    
    
    components - str,
                 path and name to the components catalogue .fits file
    
    Returns
    -------
    
    maskedComp_array - 2D array,
                       A masked array with just the source to used for the
                       ridgeline process
    """

    source_name = source[0]
    n_comp = float(source[3])  ## The number of components
    xra = float(source[4])
    ydec = float(source[5]) 
    #fits_name = "Catalogue/fits/" + source_name + '.fits'
    #numpy_name = 'Catalogue/rms4/' + source_name + '.npy'
    #maskedComp_array = flux_array
    #hdu = fits.open(fits_name)

#Create the appropriate arrays and counter to fill them in

    try: 
        n_comp=int(float(n_comp))  ## This might not be necessary so could take out at a later date.
    except ValueError:
        n_comp=1
    if int(float(n_comp)) == 0:
        n_comp = 1
        
    ell_ra=np.zeros(n_comp)
    ell_dec=np.zeros(n_comp)
    ell_maj=np.zeros(n_comp)
    ell_min=np.zeros(n_comp)
    ell_pa=np.zeros(n_comp)
    ell_flux=np.zeros(n_comp)
    compcounter = 0
    excludeComp = 0
        
    regwrite=open(RLF.tmpdir+'temp2.reg','w')  ## Creates the file for storing the infor about excluded regions
    regwrite.write('# Region file format: DS9 version 4.1\n'+'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'+'fk5\n')
        
    for row in CompTable:
        source2 = row[0]
        comp_ra = float(row[2])
        comp_dec = float(row[3])
        comp_flux = float(row[4])
        comp_maj = float(row[5])
        comp_min = float(row[6])
        comp_pa = float(row[7])
        
        if (source_name == source2):
            if n_comp>=1:
                ell_ra[compcounter]=comp_ra
                ell_dec[compcounter]=comp_dec
                ell_flux[compcounter]=comp_flux
                ell_maj[compcounter]=comp_maj
                ell_min[compcounter]=comp_min
                ell_pa[compcounter]=comp_pa
                compcounter+=1
            else:
                ell_ra=comp_ra
                ell_dec=comp_dec
                ell_flux=comp_flux
                ell_maj=comp_maj
                ell_min=comp_min
                ell_pa=comp_pa
#Defining ellipse regions to mask out unrelated components.
        else:
            try:
                compSep=1.0
                if (abs(float(comp_ra)-float(xra))<=0.4 and abs(float(comp_dec)-float(ydec))<=0.4):
                    coordsOpt=SkyCoord(ra=float(xra)*u.degree, dec=float(ydec)*u.degree, frame='fk5')
                    coordsComp=SkyCoord(ra=float(comp_ra)*u.degree, dec=float(comp_dec)*u.degree, frame='fk5')
                    compSep=coordsOpt.separation(coordsComp)/u.degree
                    if compSep<=0.5:
                        regwrite.write('ellipse('+str(comp_ra)+','+str(comp_dec)+','+str(comp_maj)+'",'+str(comp_min)+'",'+str((float(comp_pa)+90.0))+')\n')
                        if excludeComp!=1:
                            excludeComp=1
            except ValueError:
                pass
    regwrite.close()

#[OTHER CODE TO TRANSFORM COORDINATES, ETC]
#First we mask out unrelated components, ( Beas code )

    if excludeComp==1:
        exRegion=pyregion.open(RLF.tmpdir+'temp2.reg').as_imagecoord(hdu[0].header)
        exMask=exRegion.get_mask(hdu=hdu[0])
        #print(maskedComp_array.mask)
        #maskedComp_array[exMask]=np.ma.masked
        maskedComp_mask = np.ma.mask_or(flux_array.mask, exMask)
        maskedComp_array = np.ma.array(flux_array, mask = maskedComp_mask)
            
    return maskedComp_array

#############################################

def GetProblematicSources(problem_filename):

    """ 
    Returns a list of source names for the problematics sources obtained
    from the problematic source list.
    
    Parameters
    ----------
    
    problem_filename - str,
                       The filename of the problematic sources obtained 
                       when having run the full cycle of ridgeline code.
    
    Returns
    -------
    
    problem_sources - a list of source names    
    
    """
    
    problem_list = open(problem_filename, 'r').readlines()
    problem_sources=[]
    for psource in problem_list:
        psource_info = psource.strip()
        pinfo_cols = psource_info.split()
        problem_sources.append(pinfo_cols[0].strip("''"))

    return problem_sources

############################################# 

def GetRidgeCone(area_fluxes, maxima_array, ridge_point, ridge_phi, chain_mask, dphi):

    """
    Returns a new search cone along the ridge with an end point
    check cone.
    WARNING: FUNCTION UNDER DEVELOPMENT.
    Issues to be fixed - finding the source end point to avoid swirly
    ridges. Directly connected to function GetRidgePoint issue.

    Parameters
    -----------
    area_fluxes - 2D array,
                  raw image array convolved with a centre-heavy cross
                  kernel

    maxima_array - 2D array, shape(area_fluxes.shape)
                   array of boolean values with pixels containing
                   a maximum having value = True

    ridge_point - 1D array, shape(2, )
                  previous point along the ridge line, from which the
                  next point is sought

    ridge_phi - float,
                angular direction of the previous ridge point

    chain_mask - masked array mask,
                 already existing mask, which covers an R pix radius
                 around previous ridge points. In place to avoid the
                 ridge line wrapping onto itself

    dphi - float,
           1/2 of the value of the cone opening angle

    Returns
    -----------

    next_cone - 2D masked array,
                area_fluxes array masked to preserve the new
                search cone only

    check_cone - 2D masked array, (NEEDS DEVELOPMENT)
                 area_fluxes array masked to preserve the stop finding
                 cone.

    Notes
    -----------
    UNDER DEVELOPMENT.
    DEVELOPER NOTES: fix the endpoint finding criteria. (possibly
    the check cone is to blame)

    """

    r, phi = PolarCoordinates(area_fluxes, ridge_point)
    maxima_mask = np.ma.mask_or(~maxima_array, chain_mask)
    masked_r = np.ma.masked_array(r, mask=maxima_mask)

    if np.ma.count(masked_r) == 0:
        next_phirange = AngleRange(ridge_phi, np.radians(dphi),\
                                   np.radians(dphi))
        next_cone = GetCone(area_fluxes, ridge_point, next_phirange)

    else:
        maxima_phi = phi[np.unravel_index(np.ma.argmin(masked_r), \
                                          masked_r.shape)]
        quad1, quad2 = CheckQuadrant(maxima_phi), CheckQuadrant(ridge_phi)

        if (quad1 == 3 and quad2 == 2) or (quad1 == 2 and quad2 == 3):
            next_phi = PiRange((TwoPiRange(maxima_phi) + TwoPiRange(ridge_phi))/2.)
            next_phirange = AngleRange(next_phi, np.radians(dphi), np.radians(dphi))
            next_cone = GetCone(area_fluxes, ridge_point, next_phirange)

        else:
            next_phi = (maxima_phi + ridge_phi)/2.
            next_phirange = AngleRange(next_phi, np.radians(dphi),\
                                       np.radians(dphi))
            next_cone = GetCone(area_fluxes, ridge_point, next_phirange)

    check_range = AngleRange(ridge_phi, atan2(1.,6.), atan2(1.,6.))
    check_cone = GetCone(area_fluxes, ridge_point, check_range)

    return next_cone, check_cone

#############################################

def GetRidgePoint(area_fluxes, ridge_point, cone, check_cone, chain_mask, R, lmsize, Rtot, Rmax):

    """
    Returns a new point along the ridge with its angular direction
    and updated ridge chain mask (R pix radius mask around previous
    ridge points found in given direction).  This function will jump
    over gaps if the total length of the ridgeline drawn so far is 
    less that one third that of Rmax. It outputs the length of each
    ridge step taken to allow for summation in FindRidges.
    WARNING: FUNCTION UNDER DEVELOPMENT.
    Issues to be fixed - finding the source end point to avoid swirly
    ridges

    Parameters
    -----------
    area_fluxes - 2D array,
                  raw image array convolved with a centre-heavy cross
                  kernel

    ridge_point - 1D array, shape(2, )
                  previous point along the ridge line, from which the
                  next point is sought

    cone - 2D masked array,
           area_fluxes array masked to preserve the search cone only

    check_cone - 2D masked array, (NEEDS DEVELOPMENT)
                 area_fluxes array masked to preserve the stop finding
                 cone.

    chain_mask - masked array mask,
                 already existing mask, which covers an R pix radius
                 around previous ridge points. In place to avoid the
                 ridge line wrapping onto itself

    R - float,
        ridge finding step, given in pixels
        
    Rtot - float,
           current running total length of the ridgeline so far
           
    Rmax - float,
           maximum length the ridgeline can be drawn to
    
    Constants
    ---------
    
    JLim - float,
           the mulitplier of maximum ridgeline length that the line
           can jump under to find another location to draw from
    
    MaxJ - float,
           the multiplier of the source size which determines the maximum
           length the ridgeline can jump to
    
    JMax - float,
           multiplier of the ridgeline length which determines the maximum
           length the ridgeline can jump to
    

    Returns
    -----------

    new_point - 1D array, shape(2, )
                new point along the ridge line

    new_phi - float,
              angular direction of the new ridge point, found by the
              weighted average of polar coordinate angles at the given
              search point

    chain_mask - masked array mask,
                 updated input mask
                 
    new_R - float,
            the size of the ridge step

    Notes
    -----------
    UNDER DEVELOPMENT.
    DEVELOPER NOTES: fix the endpoint finding criteria

    """

    r, phi = PolarCoordinates(area_fluxes, ridge_point)
    cone_mask = cone.mask
    #check_mask  = check_cone.mask
    r_mask = np.ma.mask_or(np.ma.masked_outside(r, 0, R).mask, chain_mask)
    #r_check_mask = np.ma.mask_or(np.ma.masked_outside(r, 0, 2.).mask, chain_mask)
    slice_mask = np.ma.mask_or(r_mask, cone_mask)
    #check_slice_mask = np.ma.mask_or(r_check_mask, check_mask)
    cone_slice = np.ma.masked_array(area_fluxes, mask=slice_mask, copy=True)
    #check_slice = np.ma.masked_array(area_fluxes, mask=check_slice_mask, copy=True)

    #Rmax = np.ma.max(r)
    #Rmax = 0.75 * float(lmsize)
    safety_stop = False
    new_R = 0    
    
    if Rtot <= RLC.Jlim * Rmax:
        while np.ma.count(cone_slice) == 0:
            R += 1
            r_mask = np.ma.mask_or(np.ma.masked_outside(r, 0, R).mask, chain_mask)
            slice_mask = np.ma.mask_or(r_mask, cone_mask)
            cone_slice = np.ma.masked_array(area_fluxes, mask=slice_mask, copy = True)
            new_R = R
        
            if R > RLC.JMax * Rmax - Rtot: #RLC.MaxJ * float(lmsize)
                safety_stop = True
                new_phi = np.nan
                new_R = R
                break    
            #print('R in jump loop = ' + str(new_R))
    
    
        if safety_stop == True:

            new_point = np.full(2, np.nan)
            if RLC.debug == True:
                print('Ridgeline completed')

        else:

            cone_phislice = np.ma.masked_array(phi, mask=slice_mask, copy=True)
            quad1, quad2 = CheckQuadrant(np.ma.min(cone_phislice)), CheckQuadrant(np.ma.max(cone_phislice))
            diff = np.ma.max(cone_phislice) - np.ma.min(cone_phislice)
    
            if 3 <= quad1 <= 4 and 1 <= quad2 <= 2 and diff > pi:    
            
                for idx, phi in np.ndenumerate(cone_phislice):
                    if cone_phislice[idx] is np.ma.masked:
                        continue

                    else:
                        phi = TwoPiRange(phi)
                        cone_phislice[idx] = phi

            cone_flux = np.ma.sum(cone_slice)
            weighted_phi = np.sum(np.multiply(cone_phislice, cone_slice))
            average_phi = PiRange(weighted_phi/cone_flux)
            new_point = np.array([R*cos(average_phi), R*sin(average_phi)]) + ridge_point
            new_phi = average_phi
            new_R = R
            #print('R under Rmax value = ' + str(new_R))

        chain_mask = np.ma.mask_or(np.ma.masked_inside(r, 0, R).mask, chain_mask)


    else:
        if np.ma.count(cone_slice) == 0:

            new_point = np.full(2, np.nan)
            new_phi = np.nan
            if RLC.debug == True:
                print('Ridgeline completed')

        else:

            cone_phislice = np.ma.masked_array(phi, mask=slice_mask, copy=True)
            quad1, quad2 = CheckQuadrant(np.ma.min(cone_phislice)), CheckQuadrant(np.ma.max(cone_phislice))
            diff = np.ma.max(cone_phislice) - np.ma.min(cone_phislice)

            if 3 <= quad1 <= 4 and 1 <= quad2 <= 2 and diff > pi:    
            
                for idx, phi in np.ndenumerate(cone_phislice):
                    if cone_phislice[idx] is np.ma.masked:
                        continue

                    else:
                        phi = TwoPiRange(phi)
                        cone_phislice[idx] = phi

            cone_flux = np.ma.sum(cone_slice)
            weighted_phi = np.sum(np.multiply(cone_phislice, cone_slice))
            average_phi = PiRange(weighted_phi/cone_flux)
            new_point = np.array([R*cos(average_phi), R*sin(average_phi)]) + ridge_point
            new_phi = average_phi
            new_R = R
            #print('R over Rmax value = ' + str(new_R))

        chain_mask = np.ma.mask_or(np.ma.masked_inside(r, 0, R).mask, chain_mask)

    return new_point, new_phi, chain_mask, new_R

#############################################

def InitAngleRanges(theta1, theta2, dphi):

    """
    Returns initial search directions along with their limiting
    phi-coordinate values in the form [lower bound, upper bound]
    in the [-pi,pi] range. Takes account of search directions
    separated by less than 2*dphi, ajusting the bounds to prevent
    first-direction information from entering the second-direction
    calculations.

    Parameters
    -----------

    theta1 - float,
             first search direction in range [-pi,pi], given in radians

    theta2 - float,
             second search direction in range [-pi,pi], given in radians

    dphi - float,
           1/2 of the phi-coordinate range of the initial search cone

    Returns
    -----------

    larger - float,
             more positive of the two initial angles, given in radians

    smaller - float,
              more negative of the two initial angles, given in radians

    range_l - 1D array, shape(2,)
              angular range limiting values for the larger angle in the
              range [-pi,pi]. Ordered as [angle-dphi, angle+dphi]

    range_s - 1D array, shape(2,)
              angular range limiting values for the smaller angle in the
              range [-pi,pi]. Ordered as [angle-dphi, angle+dphi]

    Notes
    -----------
    range_l and range_s can both have more positive angle values as the
    lower limits/ negative angle values as the upper limits, when they
    exceed [-pi, pi] range.
    If theta1 and theta2 are np.nan, the ranges returned are np.nan as
    well, to allow for elimination of problematic sources.

    """

    if np.isnan(theta1) and np.isnan(theta2):
        larger, smaller = np.nan, np.nan
        range_l = range_s = np.array([np.nan,np.nan])
        print('No initial directions were found. ' \
                'Further source analysis is aborted.')

    else:

        larger = max(theta1, theta2)
        smaller = min(theta1, theta2)
        range_l, range_s = AngleRange(larger, dphi, dphi), \
                AngleRange(smaller, dphi, dphi)
        upper_lquad, lower_squad = CheckQuadrant(range_l[1]),\
                CheckQuadrant(range_s[0])

        if (upper_lquad == 2 and lower_squad == 2) or \
                (upper_lquad == 3 and lower_squad == 3) or \
                (upper_lquad == 3 and lower_squad == 2):

            if range_l[1] > range_s[0]:
                middle = PiRange((TwoPiRange(larger) + \
                                  TwoPiRange(smaller))/2.)
                range_l = np.array([range_l[0], middle])
                range_s = np.array([middle, range_s[1]])

        else:
            if range_s[1] > range_l[0]:
                middle = (larger+smaller)/2.
                range_l = np.array([middle, range_l[1]])
                range_s = np.array([range_s[0], middle])

    return larger, smaller, range_l, range_s

#############################################

def InitCones(area_fluxes, init_point, dphi, lmsize, CompTable, \
              n_comp, source_name, hdu):

    """
    Returns initial position, search directions and search
    cones ready for ridge detection.

    Parameters
    -----------

    area_fluxes - 2D array,
                  raw image array convolved with a centre-heavy cross
                  kernel

    init_point - 1D array, shape(2,)
                 starting point position for ridge detection

    dphi - float,
           1/2 of the phi-coordinate range of the initial search cone

    lmsize - str,
             the size of the source in pixels

    Returns
    -----------

    larger - float,
             more positive of the two initial search angle
             directions, given in radians

    smaller - float,
              more negative of the two initial search angle
              directions, given in radians

    cone_l - 2D array, shape(area_fluxes.shape)
             area_fluxes array masked such that unmasked values
             correspond to the first initial search cone in
             ridge detection

    cone_s - 2D array, shape(area_fluxes.shape)
             area_fluxes array masked such that unmasked values
             correspond to the second initial search cone in
             ridge detection

    Notes
    -----------
    Returns an array full of NaN if there is no initial direction
    found or the optical ID is outside of the image boundaries.

    """
    Error = 'N/A'
    maxima, maxima_array, eroded = ErodedMaxima(area_fluxes)
    #init_point = InitPoint(area_fluxes, optical_pos, maxima, lmsize)
    eroded_mask = np.ma.masked_equal(eroded, 0).mask

    # check whether the initial point exists within
    # the region of interest

    y, x = np.mgrid[slice((0),(area_fluxes.shape[0]),1),\
                    slice((0),(area_fluxes.shape[1]),1)]
    y = np.ma.masked_array(y, \
                           mask=np.ma.masked_invalid(area_fluxes).mask)
    x = np.ma.masked_array(x, \
                           mask=np.ma.masked_invalid(area_fluxes).mask)
    xmin = np.ma.min(x)
    xmax = np.ma.max(x)
    ymin = np.ma.min(y)
    ymax = np.ma.max(y)   

    if (xmin < init_point[0] < xmax) and (ymin < init_point[1] < ymax):

        theta1, theta2 = InitDirections(area_fluxes, \
                                        init_point, maxima_array)
        
        if np.isnan(theta1) and np.isnan(theta2):
            theta1, theta2 = MaskedCentre(init_point, area_fluxes)
            Error = 'No Initial Directions'
            
            q = 0  ##  Counter for iterations on initial point finder
            previp = np.full_like(init_point, init_point)
            new_iparray = np.copy(area_fluxes)
            while (np.isnan(theta1) and np.isnan(theta2)) and q < RLC.ipit:
                theta1, theta2, new_ip = RetryDirections(new_iparray, \
                                                    init_point, CompTable, n_comp, source_name, hdu)
                previp = np.vstack((previp, new_ip))
                init_point = new_ip
                if RLC.debug == True:
                    print('Retry Cones Run ' + str(q))
                new_iparray[(previp[:]).astype(int)] = np.nan
                q += 1
                
        larger, smaller, range_l, range_s = \
                        InitAngleRanges(theta1, theta2, dphi)

        if np.isnan(larger) and np.isnan(smaller):

            cone_l = cone_s = np.array([np.nan, np.nan])

        else:
            r, phi = PolarCoordinates(area_fluxes, init_point)

            if 3 <= CheckQuadrant(range_l[1]) <= 4 \
                and 1 <= CheckQuadrant(range_l[0]) <= 2:

                phil = np.ma.masked_inside(phi, range_l[1], range_l[0])
                mask_l = np.ma.mask_or(eroded_mask, phil.mask)
                cone_l = np.ma.masked_array(area_fluxes, mask=mask_l)

            else:
                phil = np.ma.masked_outside(phi, range_l[1], range_l[0])
                mask_l = np.ma.mask_or(eroded_mask, phil.mask)
                cone_l = np.ma.masked_array(area_fluxes, mask=mask_l)

            if 3 <= CheckQuadrant(range_s[1]) <= 4 \
                and 1 <= CheckQuadrant(range_s[0]) <= 2:

                phis = np.ma.masked_inside(phi, range_s[1], range_s[0])
                mask_s = np.ma.mask_or(eroded_mask, phis.mask)
                cone_s = np.ma.masked_array(area_fluxes, mask=mask_s)

            else:
                phis = np.ma.masked_outside(phi, range_s[1], range_s[0])
                mask_s = np.ma.mask_or(eroded_mask, phis.mask)
                cone_s = np.ma.masked_array(area_fluxes, mask=mask_s)

    else:
        print('Optical ID position is outside the region of interest. '\
              'Further source analysis is aborted')
        Error = 'Initial_ID_Outside_Region. Broken_Cutout'

        larger, smaller = np.nan, np.nan
        cone_l = cone_s = np.array([np.nan,np.nan])

    return larger, smaller, cone_l, cone_s, init_point, Error

#############################################

def InitDirections(area_fluxes, init_point, maxima_array):

    """
    Returns initial directions for ridge finding as angles
    in the [-pi,pi] range. Angular coordinates are calculated
    by creating a polar coordinate grid centred on the starting
    point position

    Parameters
    -----------

    area_fluxes - 2D array,
                  raw image array convolved with a centre-heavy cross
                  kernel

    init_point - 1D array, shape(2,)
                 starting point position for ridge detection

    maxima_array - 2D array, shape(area_fluxes.shape)
                   array of boolean values with pixels containing a maximum
                   having value = True

    Returns
    -----------

    theta1 - float,
             first initial direction for ridge finding, given in radians

    theta2 - float,
             second initial direction for ridge finding, given in radians

    Notes
    -----------

    Returns theta1, theta2 = np.nan, when only one maximum is detected
    and identified as the initial point. This allows for excluding
    problematic sources from the analysis in a sequence.

    """

    r0, phi0 = PolarCoordinates(area_fluxes, init_point)

    r0 = np.ma.masked_array(r0, ~maxima_array)
    phi0 = np.ma.masked_array(phi0, ~maxima_array)

    init_pix = np.floor(init_point).astype('int')
    r0[init_pix[1], init_pix[0]] = np.ma.masked
    phi0[init_pix[1], init_pix[0]] = np.ma.masked

    if np.ma.count(r0) == 0:
        theta1 = np.nan
        theta2 = np.nan

    else:
        idx1 = np.unravel_index(np.argmin(r0), r0.shape)
        theta1 = phi0[idx1]
        phirange = AngleRange(theta1, np.radians(60), np.radians(60))

        if 3 <= CheckQuadrant(phirange[1]) <= 4 and \
                1 <= CheckQuadrant(phirange[0]) <= 2:
            phi0 = np.ma.masked_outside(phi0, phirange[1], phirange[0])

        else:
            phi0 = np.ma.masked_inside(phi0, phirange[0], phirange[1])

        r0 = np.ma.masked_array(r0, mask=phi0.mask)

        if np.ma.count(r0) == 0:
            if theta1 >= 0.:
                theta2 = theta1 - pi
            else:
                theta2 = theta1 + pi

        else:
            idx2 = np.unravel_index(np.argmin(r0), r0.shape)
            theta2 = phi0[idx2]

    return theta1, theta2

#############################################

def InitPoint(area_fluxes):   

    """Returns position of the starting point for ridge detection based
    on the highest flux in the associated components in the LOFAR
    catalogue. As the area_fluxes image is already masked for the
    associated components, we just need to find the brightest point,
    hence this version of the code does not take any catalogue
    entries as input.

    Returns
    -----------

    init_point - 1D array, shape(2,)
                 starting point position for ridge detection

    Notes
    -----------
    Radius of possible offset is given by the beam size = 5 pixels IF NEEDED

    """
    maxflux=np.unravel_index(np.nanargmax(area_fluxes),area_fluxes.shape)
    return np.array([maxflux[1],maxflux[0]])
    
#############################################

def MaskedCentre(init_point, area_fluxes):
    
    """
    Masks an area around the initial starting point, when two starting
    directions cannot be found inorder to try and find two starting
    directions without a bright central concentration of flux.
    
    Parameters
    ----------
    
    init_point - 1D array, shape(2.)
                 centre pixel position
                   
    area_fluxes - array,
                  the array that corresponds to the source
    
    Constants
    ---------
    
    Rad - float,
          the size of the radius of the circle used to mask the initial
          point

    Returns
    -------
    
    NewAngle1 - float,
                The next attempt at obtaining initial directions for the first 
                ridgeline, in radians.
    
    NewAngle2 - float,
                The next attempt at obtaining initial directions for the second
                ridgeline, in radians.
    
    """
    
    masking_area = np.copy(area_fluxes)
    
    radius = RLC.Rad ## Set the radius of the masking circle
    mysize, mxsize = masking_area.shape ## Find the size of the array
    xip = init_point[0]  ## Find the x coordinate of the initial point
    yip = init_point[1]  ## Find the y coordinate of the initial point

    xcore, ycore = np.ogrid[-yip : mysize - yip, -xip : mxsize - xip]  ## Form a grid of differences between inital points and size
    circle_mask = xcore * xcore + ycore * ycore <= radius * radius  ## Draw a circle of all elements inside the region

    masking_area[circle_mask] = np.nan ## Set the values in the circle mask to nan

    maskedmaxima, maskedmaxima_array, maskederoded = ErodedMaxima(masking_area) ## Finding the maximas on the eroded array without the centre point
    
    maskedinit_point = np.empty(2)

    # Find the maximum flux point in the cutout to use as the starting point
    maskedmaxflux = np.unravel_index(np.argmax(np.ma.masked_invalid(masking_area)), masking_area.shape)
    maskedinit_point[0] = maskedmaxflux[1]
    maskedinit_point[1] = maskedmaxflux[0]

    #maskedinit_point = InitPoint(masking_area, init_point, maskedmaxima, lmsize) ## Finding the starting point with out the centre    
    
    # check whether the initial point exists within
    # the region of interest

    y, x = np.mgrid[slice((0),(masking_area.shape[0]),1),\
                    slice((0),(masking_area.shape[1]),1)]
    y = np.ma.masked_array(y, \
                           mask=np.ma.masked_invalid(masking_area).mask)
    x = np.ma.masked_array(x, \
                           mask=np.ma.masked_invalid(masking_area).mask)
    xmin = np.ma.min(x)
    xmax = np.ma.max(x)
    ymin = np.ma.min(y)
    ymax = np.ma.max(y)

    if (xmin < maskedinit_point[0] < xmax) and (ymin < maskedinit_point[1] < ymax):
        
        NewAngle1, NewAngle2 = InitDirections(masking_area, maskedinit_point, maskedmaxima_array)
        
    else:
        
        NewAngle1, NewAngle2 = np.nan, np.nan
        
    return NewAngle1, NewAngle2


#############################################

def PiRange(angle):

    """
    Ensures the angle belongs to [-pi,pi] range

    Parameters
    -----------

    angle - float,
            angle of interest, given in radians

    Returns
    -----------

    angle_out - float,
                input angle in the range [-pi,pi]
    """

    if angle < -pi:
        angle_out = angle + 2*pi

    elif angle > pi:
        angle_out = angle - 2*pi

    else:
        angle_out = angle

    return angle_out

#############################################

def PolarCoordinates(a, pos):

    """
    Returns a grid of r and phi polar coordinate values
    corresponding to each pixel, with the origin position
    provided as a parameter

    Parameters
    -----------

    a - 2D array,
        array of pixels to be mapped onto polar coordinates

    pos - 1D array, shape(2, )
          array containing central point position given as pixel
          values [x_position, y_position]

    Returns
    -----------

    r - 2D array, shape(a.shape)
        array of r coordinate values corresponding to each pixel
        in provided array a

    phi - 2D array, shape(a.shape)
          array of phi coordinate values corresponding to each
          pixel in provided array a
    """

    X = pos[0]
    Y = pos[1]
    r = np.zeros_like(a)
    phi = np.zeros_like(a)

    for index, val in np.ndenumerate(a):

        # radius value meshgrid
        r_x = (index[1]+0.5) - X
        r_y = (index[0]+0.5) - Y
        r[index[0],index[1]] = sqrt(r_x**2 + r_y**2)

        # phi value grid for [-pi,pi] range
        theta = round(atan2(r_y, r_x), 3)
        theta = PiRange(theta)
        phi[index[0], index[1]] = theta

    return r, phi

#############################################

def RetryDirections(new_iparray, init_point, CompTable, n_comp, source_name, hdu):
    
    """
    """
    
    #new_ip = np.empty(2)

    new_iparray[int(init_point[1]),int(init_point[0])] = np.nan
    #print(new_iparray)
    maskedmaxima2, maskedmaxima_array2, maskederoded2 = ErodedMaxima(new_iparray)
       
    #newmax = np.unravel_index(np.argmax(np.ma.masked_invalid(new_iparray[int(init_point[1]-1):int(init_point[1]+2), \
                                                                   #int(init_point[0]-1):int(init_point[0]+2)])), new_iparray.shape)
    new_ip = InitPoint(new_iparray)
    #print(newmax)
    #new_ip[0] = newmax[1]
    #new_ip[1] = newmax[0]
       
    NewAngle3, NewAngle4 = InitDirections(new_iparray, new_ip, maskedmaxima_array2)
    
    if np.isnan(NewAngle3) and np.isnan(NewAngle4):
        NewAngle3, NewAngle4 = MaskedCentre(new_ip, new_iparray)
    
    return NewAngle3, NewAngle4, new_ip

#############################################

def TotalFluxSelector(catalogue1, ffo):

    """
    Saves a .txt file containing source names and their LOFAR catalogue
    positions, number of components, size of source in 
    arcsecs, total flux and Host ID location in separate rows.  This is
    based on the total flux selection criterion from the
    eyeclass_4rms_v1.2_morph_agn_length_err.fits file.
    Also allows for the binning of this data into different bins based on
    size, flux and number of LOFAR Galaxy Zoo associated components.

    Parameters
    -----------

    catalogue1 - str,
                 path and name to the catalogue .fits file


    ffo: a Flood object from sizeflux
    
    Constants
    ---------
    
    LGZA - int,
           if active, number of LOFAR Galaxy Zoo components associated with
           the source.  Needs to the dot to do the filtering.

    Notes
    -----------
    The structure of .txt file is based on the sarah_sddsmatches.txt.

    """
    print("Total flux selector opening: ",catalogue1)
    #print("Components table is ",CompTable)
    tbdata = Table.read(catalogue1)  ## Open the radio Catalogue
    tbdata['Size'] = np.where(~np.isnan(tbdata['LGZ_Size']),tbdata['LGZ_Size']/3600.0,2*tbdata['Maj']*2.0/3600.0)  ## Find the source size in deg -- gives estimate to filter exclusion catalogue only
    print("First source is ",tbdata[0][RLF.SSN])
    sizes = []
    for r in tbdata:
        source_name = r[RLF.SSN]
        Lra = r[RLF.SRA]
        Ldec = r[RLF.SDEC]
        n_comp = r[RLF.SASS]
        isize=r['Size']
        if n_comp == 0:
            n_comp = 1
        
        hdu = DefineHDU(source_name)
        flux_array = GetFluxArray(source_name)
        print("Max of flux array is ",np.nanmax(flux_array))
        hdu[0].data=flux_array # otherwise not passed to mask!
        cinc,cexc=ffo.select(source_name,Lra,Ldec,isize,verbose=False)

        try:
            _,FMArray = ffo.mask(source_name, cinc, cexc, hdu, None, verbose=False)
            #FMArray = FloodMask(source_name, Lra, Ldec, n_comp, CompTable, hdu, flux_array)
            print("Max of FMArray is ",np.nanmax(FMArray))
        except:
            print(source_name, 'Flood Fill and Mask Error. Source Not Suitable')
            size = np.nan
            sizes.append(size)
            pass
        else:
            size = FindSizeOfSource(FMArray)
            print(source_name, size)
            sizes.append(size)
        print("Final size and flux: ",size,r['Total_flux'])

    tbdata['Size']=sizes
    print(tbdata)

    # Note that although a Total_flux criterion is retained flux is
    # *not measured* in the above and by definition will always be >
    # 0. So this boils down to whether the preliminary mask can
    # measure a size. We could measure total flux but don't gain much...
    
    sub5 = tbdata['Size'] > 0 ## S2  THE ONE I WANT Need to have values to remove any NaN
    sub9 = tbdata['Total_flux'] > 0 ## F4  THE ONE I WANT Need to have values to remove any NaN

    selection1 = sub5 & sub9
    
    ##  Filters for the flux, size and components
    #sub10 = table['Assoc'] == RLC.LGZA
    #selection2 = [n & o for n, o in zip(selection1, sub10)]  
    
    output = tbdata[selection1]  ##  THE ONE I WANT

    return output

#############################################

def TrialSeries(available_sources, components, R, dphi, ffo):

    """
    Performs ridge finding and output plotting for a list of
    available sources with given ridge finding step and search
    radius.
    Outputs source images with overplotted ridge lines as well as
    a .txt files with ridge point positions, angular directions
    and length of ridgeline at each point.
    In case of problematic sources outputs the images with LOFAR
    catalogue positions and initial positions overplotted as well as
    a list of probelmatic sources along with the identified issue in
    a .txt file.
    Also allows for a sample of the sources to be taken.
    **NOTE - Matplotlib is currently throwing up a warning about the
    colour fill on pcolormesh.  Hence the simple filter warning ignore
    in the notebook.**

    Parameters
    -----------

    available_sources - 2D array, shape(n, 3), dtype='str'
                        array containing n source names and optical ID
                        locations in pixel values for all files available
                        in a given directory
    
    components - str,
                 path and name to the components catalogue .fits file
    
    R - float,
        ridge finding step, given in pixels

    dphi - float,
           1/2 of the value of the cone opening angle
    
    ffo - Flood object
    
    Constants
    ---------
    
    Start - int,
            the index to start the selection of sources from (default = 0)
    
    Step - int,
           the fraction step size of the function. 1 will run everything,
           2 half etc
    
    ImSize - float,
             the fraction of the source size that the final image is cut
             down to. 1 will match the size of the source from the centre
             of the cutout

    Notes
    -----------
    It is the first trial series function, which does not include automated
    classification into bent and straight sources.

    Folder structure required:
    - dataset folder containing: catalogue file + .csv file or ready .txt file
    with catalgoue matches, this toolkit
        - SAMPLE
            -4rms
            -4rms_cutouts
            -fits
            -fits_cutouts
            -ridges
            -problematic

    """

    palette = plt.cm.cividis
    palette = copy.copy(plt.cm.get_cmap("cividis"))
    palette.set_bad('k',0.0)

    problem_names = np.array([str('#Source_Name'), str('Problem_Type')])
    source_counter = 1
    Error = 'N/A'

    # To take every source
    #for source in available_sources:
    
    # To get every other source (ie 50% subset)
    #for source in available_sources[0::2]:
    
    # To get every third source (ie 33.3% subset)
    #for source in available_sources[0::3]:
    
    # To get every fourth source (ie 25% subset)
    #for source in available_sources[0::4]:
    
    # To get a sample of a certain size (ie [0::n] gives every nth source)
    for source in available_sources[RLC.Start::RLC.Step]:
        
        source_name = source[RLF.SSN]
        n_comp = source[RLF.SASS]
        Lra = source[RLF.SRA]
        Ldec = source[RLF.SDEC]
        lmsize = source['Size']
        flux_array = GetCutoutArray(source_name)
        hdu = DefineCutoutHDU(source_name)
        #hdu = DefineHDU(source_name)
        #flux_array = GetFluxArray(source_name)
        hdu[0].data=flux_array # otherwise not passed to mask!
        
        optical_pos = (float(lmsize), float(lmsize))
        
        if n_comp == 0:
            n_comp = 1
        
        print('-------------------------------------------')
        print('Source Number %s' %source_counter)
        print('Source %s under ridgeline analysis' %source_name)
        
        #area_fluxes1 = AreaFluxes(flux_array)      ## Convolution  
        #rms_array = GetRMSArray(source, hdu, flux_array)
        if RLC.debug == True:
            print(str(source_name)+' Flood Fill and Masking')
        try:
            cinc,cexc=ffo.select(source_name,Lra,Ldec,lmsize/3600.0)
            _,FandM_array=ffo.mask(source_name,cinc,cexc,hdu,None,verbose=False)
        except IndexError:
            # This try/except is not needed now but preserved to match
            # old code structure
            print('Flood fill error')
            raise
        else:
            #print(str(source_name)+' Masking')
            #maskedComp_array = GetMaskedComp(hdu, source, components, flooded_array, CompTable)  ## Masking     
            if RLC.debug == True:
                print(str(source_name)+' Convolution')
            area_fluxes = AreaFluxes(FandM_array)  ## Convolution
                    
            if RLC.debug == True:
                print(str(source_name)+' Initial Point')
            
            init_point = InitPoint(area_fluxes)
            
            if RLC.debug: print('Init point is',init_point)

            # find both ridges and angular directions at each step
            if RLC.debug == True:
                print(str(source_name)+' Ridgeline')
            try:
                ridge1, phi_val1, Rlen1, ridge2, phi_val2, Rlen2, Error = FindRidges(area_fluxes, \
                                    init_point, R, dphi, lmsize, ffo.t, n_comp, source_name, hdu)
                    
            except (ValueError, UnboundLocalError):#,TypeError, 
            
                    y, x = np.mgrid[slice((0),(area_fluxes.shape[0]),1), slice((0),(area_fluxes.shape[1]),1)]
                    fig, ax = plt.subplots(figsize=(8,8))
                    fig.suptitle('Source: %s' %source_name)
                    fig.subplots_adjust(top=0.9)
                    ax.set_aspect('equal', 'datalim')
                    A = np.ma.array(area_fluxes, mask=np.ma.masked_invalid(area_fluxes).mask)
                    ax.pcolor(x, y, A, cmap=palette, vmin=np.nanmin(A), vmax=np.nanmax(A))
                    ax.scatter(float(optical_pos[0]), float(optical_pos[1]), s=130, c='m', marker='x', label='LOFAR id')
                    ax.scatter(float(init_point[0]), float(init_point[1]), s=130, c='c', marker='x', label='Initial point')
                    ax.set_xlim(x.min(), x.max()), ax.set_ylim(y.min(), y.max())
                    ax.legend()
                                           
                    problem = np.array([str(source_name),str('Initial_Directions_Error_Occurred')])
                    problem_names = np.vstack((problem_names, problem))
                    fig.savefig(RLF.Probs %source_name)
                    plt.close(fig)
                    print('Error Occurred. No Initial Directions Found. Further Source Analysis is Aborted.')
            
            except (TypeError):#UnboundLocalError,TypeError, 
            
                    y, x = np.mgrid[slice((0),(area_fluxes.shape[0]),1), slice((0),(area_fluxes.shape[1]),1)]
                    fig, ax = plt.subplots(figsize=(8,8))
                    fig.suptitle('Source: %s' %source_name)
                    fig.subplots_adjust(top=0.9)
                    ax.set_aspect('equal', 'datalim')
                    A = np.ma.array(area_fluxes, mask=np.ma.masked_invalid(area_fluxes).mask)
                    ax.pcolor(x, y, A, cmap=palette, vmin=np.nanmin(A), vmax=np.nanmax(A))
                    ax.scatter(float(optical_pos[0]), float(optical_pos[1]), s=130, c='m', marker='x', label='LOFAR id')
                    ax.scatter(float(init_point[0]), float(init_point[1]), s=130, c='c', marker='x', label='Initial point')
                    ax.set_xlim(x.min(), x.max()), ax.set_ylim(y.min(), y.max())
                    ax.legend()
                                           
                    problem = np.array([str(source_name),str('Erosion_Error_Occurred')])
                    problem_names = np.vstack((problem_names, problem))
                    fig.savefig(RLF.Probs %source_name)
                    plt.close(fig)
                    print('Erosion Error Occurred. Further Source Analysis is Aborted.')            
            
            else:
                
                    if np.all(np.isnan(ridge1)) and np.all(np.isnan(ridge2)):
    
                        y, x = np.mgrid[slice((0),(area_fluxes.shape[0]),1), slice((0),(area_fluxes.shape[1]),1)]
                        fig, ax = plt.subplots(figsize=(8,8))
                        fig.suptitle('Source: %s' %source_name)
                        fig.subplots_adjust(top=0.9)
                        ax.set_aspect('equal', 'datalim')
                        A = np.ma.array(area_fluxes, mask=np.ma.masked_invalid(area_fluxes).mask)
                        ax.pcolor(x, y, A, cmap=palette, vmin=np.nanmin(A), vmax=np.nanmax(A))
                        ax.scatter(float(optical_pos[0]), float(optical_pos[1]), s=130, c='m', marker='x', label='LOFAR id')
                        ax.scatter(float(init_point[0]), float(init_point[1]), s=130, c='c', marker='x', label='Initial point')
                        ax.set_xlim(x.min(), x.max()), ax.set_ylim(y.min(), y.max())
                        ax.legend()
                        
                        #problem = np.array([str(source_name),str('%s_Error_Occurred' %Error)])
                        #problem_names = np.vstack((problem_names, problem))                       
        
                        if np.any(np.isnan(optical_pos)):
    
                            problem = np.array([str(source_name), str('Optical_ID_issue')])
                            problem_names = np.vstack((problem_names, problem))
                            print('Optical position evaluates to NaN')
    
                        elif x.min()<= float(optical_pos[0])<=x.max() and y.min()<=float(optical_pos[1])<=y.max():
    
                            problem = np.array([str(source_name), str('%s. Initial_direction_issue' %Error)]).astype('str')
                            problem_names = np.vstack((problem_names, problem))
                            #ax.scatter(float(optical_pos[0]), float(optical_pos[1]), c='m', s=45, label='LOFAR id')
                            #ax.scatter(float(init_point[0]), float(init_point[1]), s=130, c='c', marker='x', label='Initial point')
                            ax.legend()

                        else:
    
                            problem = np.array([str(source_name), str('Optical_ID_issue')])
                            problem_names = np.vstack((problem_names, problem))
                            print('Optical id issue: can be out of the image')
    
                        fig.savefig(RLF.Probs %source_name)
                        plt.close(fig)
    
                    else:
    
                # save the output as a .txt file of the form:
                # x_coord y_coord angle_dir (3 columns, space separated)
                        file1 = np.column_stack((ridge1[:,0], ridge1[:,1], phi_val1, Rlen1))
                        file2 = np.column_stack((ridge2[:,0], ridge2[:,1], phi_val2, Rlen2))
                        np.savetxt(RLF.R1 %source_name, file1, delimiter=' ')
                        np.savetxt(RLF.R2 %source_name, file2, delimiter=' ')
    
                        y, x = np.mgrid[slice((0),(area_fluxes.shape[0]),1), slice((0),(area_fluxes.shape[1]),1)]
                        y = np.ma.masked_array(y, mask=np.ma.masked_invalid(area_fluxes).mask)
                        x = np.ma.masked_array(x, mask=np.ma.masked_invalid(area_fluxes).mask)
                        
                        xmin = np.ma.min(x)
                        xmax = np.ma.max(x)
                        ymin = np.ma.min(y)
                        ymax = np.ma.max(y)
                        
                        x_source_min = float(optical_pos[0]) - RLC.ImSize * float(lmsize)
                        x_source_max = float(optical_pos[0]) + RLC.ImSize * float(lmsize)
                        y_source_min = float(optical_pos[1]) - RLC.ImSize * float(lmsize)
                        y_source_max = float(optical_pos[1]) + RLC.ImSize * float(lmsize)
                        
                        if x_source_min < xmin:
                            xplotmin = xmin
                        else:
                            xplotmin = x_source_min
                                
                        if x_source_max < xmax:
                            xplotmax = x_source_max
                        else:
                            xplotmax = xmax
                        
                        if y_source_min < ymin:
                            yplotmin = ymin
                        else:
                            yplotmin = y_source_min
                                
                        if y_source_max < ymax:
                            yplotmax = y_source_max
                        else:
                            yplotmax = ymax
                        
                        
                        fig, ax = plt.subplots(figsize=(10,10))
                        fig.suptitle('Source: %s' %source_name)
                        fig.subplots_adjust(top=0.9)
                        ax.set_aspect('equal', 'datalim')
    
                        A = np.ma.array(area_fluxes, mask=np.ma.masked_invalid(area_fluxes).mask)
                        ax.pcolor(x, y, A, cmap=palette, vmin=np.nanmin(A), vmax=np.nanmax(A))
                        ax.plot(ridge1[:,0], ridge1[:,1], 'r-', label='ridge 1')
                        ax.plot(ridge2[:,0], ridge2[:,1], 'r-', label='ridge 2')
                        #ax.scatter(float(cat_pos[0]), float(cat_pos[1]), s=130, c='m', marker='x', label='LOFAR id')
                        #ax.scatter(float(init_point[0]), float(init_point[1]), s=130, c='c', marker='x', label='Initial point')
                        ax.legend()
                        ax.set_xlim(xplotmin, xplotmax)
                        ax.set_ylim(yplotmin, yplotmax)
    
                        fig.savefig(RLF.Rimage %(source_name, dphi))
                        plt.close(fig)

    # save the list of problematic sources in problematic_sources_list.txt
    # with source names separated from the problem type with a space.
    # the file has a header beginning with a hash # which describes the
    # content of each column

        np.savetxt(str(RLF.psl), problem_names, fmt='%s', delimiter=' ')
        source_counter += 1
        
#############################################

def TwoPiRange(angle):

    """
    Ensures the angle belongs to [0,2pi] range

    Parameters
    -----------

    angle - float,
            angle of interest, given in radians

    Returns
    -----------

    angle_out - float,
                input angle in the range [0,2pi]
    """

    if angle < 0:
        angle_out = 2*pi + angle

    else:
        angle_out = angle

    return angle_out

##########################################################################################


##  THESE FUNCTIONS ARE REDUNDANT BUT I HAVE KEPT THEM IN ANYWAY

def FindNoise(source_array):
    
    """
    Returns the mean and standard deviation of an array.
    
    Parameters
    ----------
    
    source_array - 2D array,
                   the original array for the source
    
    
    Returns
    -------
    
    source_array_mean - float,
                        the mean of the array
                   
    source_array_std - float,
                       the standard deviation of the array
    
    """
    
    flat_source_array=source_array.flatten()
    source_array_mean=np.nanmean(flat_source_array)
    source_array_std=np.nanstd(flat_source_array)
    return source_array_mean,source_array_std

#############################################


def FindNoiseArea(source, hdu3):
        
    """
    Finds an area and then works out the noise on that area
    
    Parameters
    ----------
    
    source - str,
             the source from available sources that the noise is
             to be found for
    
    hud3 - .fits
           the loaded .fits file for the source
    
    
    Returns
    -------
    
    mean - float,
           the mean of the array
    
    noise - float,
            the noise of the array
    
    ## Working Notes - this function is mostly un-needed and un used
    as there is not a size cutout of the whole cut out used, but I 
    have kept it in as part of the cycle ##
    
    """

    source_name = source[0]  
    ysize2,xsize2=hdu3[0].data.shape
    print('-------------------------------------------')
    print('Source %s under noise analysis' %source_name)
    print("ysize,xsize: ",ysize2,xsize2)
    #xmin2, ymin2 = [0,0]
    #xmax2, ymax2 = xsize2, ysize2
    #print(xmin2,xmax2,ymin2,ymax2)

    subim = hdu3[0].data[0:ysize2,0:xsize2]
    mean,noise = FindNoise(subim)

    print("mean, noise: ", mean, noise)
    return mean, noise

#############################################

def GetRMSArray(possible_source, hdu4, flux_array):
    
    """
    Returns an array with the noise reduction applied to it
    
    Parameters
    ----------
    
    possible_source - str,
                      the possible source from available sources
    
    hdu4 - .fits
           the loaded .fits file
    
    flux_array - 2D array
                 the array that corresponds to the source
    
    Constants
    ---------
    
    nSig - int,
           the multiple of sigma for RMS reduction
    
    Returns
    -------
    
    rmsed_array - 2D array,
                  array that has had the noise removed from it
    
    ## Working notes - this has currently been removed as it interfers 
    with the erosion code however I have kept it incase it might be 
    included again at a later date. ##
    
    """
    rmsed_array = flux_array
    rms = FindNoiseArea(possible_source, hdu4)[1]
    #print(rms)
    thres = RLC.nSig * rms
    #print(thres)
    rmsed_array[rmsed_array<thres]=np.nan
    #print(flux_array)
    #print(rmsed_array)
    return rmsed_array

#############################################
   
# Length function
# Loop over x and y values if greater than 0
# Print maximum value
def length(array):
    """
    

    Parameters
    ----------
    array : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
   
    sizy, sizx = array.shape
    ddd = np.zeros_like(array)
    for p in range(sizx):
        for q in range(sizy):
            if not np.isnan(array[q,p]):
                ddd[q,p] = maxdist(array,p,q)
    return np.max(ddd)

########################

def FindSizeOfSource(array):
    """
    

    Parameters
    ----------
    source : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    #print(source)
    #source.map(lambda x: x.strip('b').strip("''"))
    #print(source)
    #array = GetFluxArray(str(source))
    
    sourcesize = RLC.ddel*3600.0*length3(array)

    return sourcesize 
#############################################

def maxdist(array,x,y):
    """
    

    Parameters
    ----------
    a : TYPE
        DESCRIPTION.
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.

    Returns
    -------
    dmax : TYPE
        DESCRIPTION.

    """
    # Set up xv and yv, make grid, calculate dx and dy, find distances, return maximum distance
    ny,nx = array.shape
    xa = np.linspace(0, nx-1, nx)
    ya = np.linspace(0, ny-1, ny)
    xv, yv = np.meshgrid(xa, ya)
    deltax = xv-x
    deltay = yv-y
    d2dx = deltax * deltax
    d2dy = deltay * deltay
    d2dxdy = d2dx + d2dy
    ddxdy = np.sqrt(d2dxdy)
    dmax = np.max(ddxdy[~np.isnan(array)])
    
    return dmax 


#############################################

# Flood-filling and masking function (uses scikit image label)
def FloodMask(source_name, xra, ydec, n_comp, CompTable, hdu, flux_array):
    
    #source_name = source[0]
    #n_comp = source[3]  ## The number of components
    #xra = float(source[4])
    #ydec = float(source[5]) 
    #fits_name = "Catalogue/fits/" + source_name + '.fits'
    #numpy_name = 'Catalogue/rms4/' + source_name + '.npy'
    #maskedComp_array = flux_array
    #hdu = fits.open(fits_name)

    #Create the appropriate arrays and counter to fill them in
    flooded_array = flux_array.copy()
    
    print("FM: source and ncomp is ",str(source_name),str(xra),str(ydec))

    try: 
        n_comp=int(float(n_comp))  ## This might not be necessary so could take out at a later date.
    except ValueError:
        n_comp=1
    #if n_comp == float(0.0):
    #    n_comp = float(1.0)
    ell_ra=np.zeros(n_comp)
    ell_dec=np.zeros(n_comp)
    ell_maj=np.zeros(n_comp)
    ell_min=np.zeros(n_comp)
    ell_pa=np.zeros(n_comp)
    ell_flux=np.zeros(n_comp)
    compcounter = 0
    excludeComp = 0
        
    regwrite=open(RLF.tmpdir+'/temp5.reg','w')  ## Creates the file for storing the infor about excluded regions
    regwrite.write('# Region file format: DS9 version 4.1\n'+'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'+'fk5\n')
    print("Opening comp table of length ",len(CompTable))   
    for row in CompTable:
        #print("processing row")
        #print("FM: Source is ",source_name,"X")
        #print(row)
        
        source2 = row[0]
        comp_ra = float(row[2])
        comp_dec = float(row[3])
        comp_flux = float(row[4])
        comp_maj = float(row[5])
        comp_min = float(row[6])
        comp_pa = float(row[7])
        
        #source2 = row['Parent_Source']
        source2 = source2.rstrip()
        #print("FM: source2 is ",source2,"X")
        if (source_name == source2):
            #print("FM: Found a match, ncomp = ",n_comp)
            if n_comp>=1:
                ell_ra[compcounter]=comp_ra
                ell_dec[compcounter]=comp_dec
                ell_flux[compcounter]=comp_flux
                ell_maj[compcounter]=comp_maj
                ell_min[compcounter]=comp_min
                ell_pa[compcounter]=comp_pa
                compcounter+=1
            else:
                ell_ra=comp_ra
                ell_dec=comp_dec
                ell_flux=comp_flux
                ell_maj=comp_maj
                ell_min=comp_min
                ell_pa=comp_pa
            
#Defining ellipse regions to mask out unrelated components.
        else:
            try:
                compSep=1.0
                #print("ellipse props: ",comp_ra,comp_dec,comp_maj,comp_min,comp_pa)
                if (abs(float(comp_ra)-float(xra))<=0.4 and abs(float(comp_dec)-float(ydec))<=0.4):
                    coordsOpt=SkyCoord(ra=float(xra)*u.degree, dec=float(ydec)*u.degree, frame='fk5')
                    coordsComp=SkyCoord(ra=float(comp_ra)*u.degree, dec=float(comp_dec)*u.degree, frame='fk5')
                    compSep=coordsOpt.separation(coordsComp)/u.degree
                    if compSep<=0.5:
                        #print("Found a close component and writing it to temp5 - ra is",comp_ra)
                        regwrite.write('ellipse('+str(comp_ra)+','+str(comp_dec)+','+str(comp_maj)+'",'+str(comp_min)+'",'+str((float(comp_pa)+90.0))+')\n')
                        if excludeComp!=1:
                            excludeComp=1
            except ValueError:
                #print("Found value error?")
                pass
    regwrite.close()
    #print("Finished writing to temp5")
#[OTHER CODE TO TRANSFORM COORDINATES, ETC]
#First we mask out unrelated components, ( Beas code )

    #if excludeComp==1:
    #    exRegion=pyregion.open('temp2.reg').as_imagecoord(hdu[0].header)
    #    exMask=exRegion.get_mask(hdu=hdu[0])
        #print(maskedComp_array.mask)
        #maskedComp_array[exMask]=np.ma.masked
    #    maskedComp_mask = np.ma.mask_or(flux_array.mask, exMask)
    #    maskedComp_array = np.ma.array(flux_array, mask = maskedComp_mask)
    
    
    #Preliminary tasks
    #Create temporary ds9 region file to use as initial mask
    regwrite=open(RLF.tmpdir+'/temp6.reg','w')
    regwrite.write('# Region file format: DS9 version 4.1\n'+'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'+'fk5\n')
    #Define auxiliary array filled with "1", which we will use to create the masks, and auxiliary data array to fill in the ellipses 
    sizey,sizex=flux_array.shape
    one_mask=1*np.isfinite(np.zeros((sizey,sizex)))
    for n in range(0,n_comp):
        regwrite.write('ellipse('+str(ell_ra[n])+','+str(ell_dec[n])+','+str(ell_maj[n])+'",'+str(ell_min[n])+'",'+str(ell_pa[n]+90.0)+')\n')
    regwrite.close()
    #In the auxiliary data array, turn to an arbitrary non-zero value the pixels inside the ellipse
    region=pyregion.open(RLF.tmpdir+'/temp6.reg').as_imagecoord(hdu[0].header)
    mask=region.get_mask(hdu=hdu[0])
    flux_array[mask==1]=0.02
    #Masking out nearby sources (previously in main code)
    if excludeComp>0:
        region2=pyregion.open(RLF.tmpdir+'/temp5.reg').as_imagecoord(hdu[0].header)
        maskOut=region2.get_mask(hdu=hdu[0])
        #Making sure that we only mask where the cutout is defined
        exclude_overlap=maskOut+one_mask        
        flux_array[exclude_overlap==2]=np.nan
    #Transform the data array into boolean, then binary values
    data_bin=1*(np.isfinite(flux_array))
    #Label the data array: this method will assign a numerical label to every island of "1" in our binary data array; we want 8-style connectivity
    data_label=label(data_bin, connectivity=2)
    #Multiply the label array by the source regions array (post-masking of extraneous sources). This allows us to identify which labels correspond to the clusters of pixels belonging to the source
    include_overlap=data_label*mask
    #Get the list of label values, excluding zero
    multi_labels=np.unique(include_overlap[np.nonzero(include_overlap)])
    #print str(multi_labels)
    #Initialise the cumulative mask
    multi_mask=np.zeros((sizey,sizex))
    #Main loop
    for i in range (0,len(multi_labels)):
        #Because of how masked arrays work, we need to explicitly set to zero the areas of the array we want masked out... we need a temporary masked array for the intermediate step
        temp_mask=(np.ma.masked_where(data_label!=multi_labels[i],one_mask))
        temp_mask[temp_mask!=1]=0
        #As we have used a 1/0 matrix as a basis, iteratively adding the island masks together will give us the full mask we need 
        multi_mask=(multi_mask+temp_mask)
    #The output mask will only contain "1" in the areas corresponding to the islands of interest
    flooded_mask=multi_mask
    #palette=plt.cm.viridis
    #palette.set_bad('k',0.0)
    #plt.rcParams['figure.figsize']=(10.67,8.0)
    #A=np.ma.array(data_label,mask=np.isnan(data_label))
    #y,x = np.mgrid[slice((0),(sizey),1),
    #               slice((0),(sizex),1)]
    #plt.pcolor(x,y,A,cmap=palette,vmin=np.nanmin(A),vmax=np.nanmax(A))
    #plt.rc('xtick', labelsize=14)
    #plt.rc('ytick', labelsize=14)
    #plt.axis([x.min(),x.max(),y.min(),y.max()]) 
    #plt.colorbar()
    #plt.savefig('Test2.pdf', dpi=300, format='pdf')
    #plt.clf()
    
    flooded_array[flooded_mask==0]=np.nan
    
    return(flooded_array)

#############################################

