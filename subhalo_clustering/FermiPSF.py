import numpy as np
import pyfits


def PSF_130(convType = 'front'):
    '''
    Generates PSF table at 130 GeV for 4-year Fermi data in GC region.
    
    Inputs: 
        -convType: Event conversion type.  Either 'front' or 'back' 
    Return:
     -theta: list of radii
     -PSF: Probabilities normalized such that max(PSF)=1.0
    '''     
    #===========================================================================
    # Load PSF from FITS averaged over current GC fermi data
    #===========================================================================
    hdulist = 0
    if convType == 'front':
        hdulist = pyfits.open('./psfFront.txt', mode='update')
    elif convType == 'back': 
        hdulist = pyfits.open('./psfBack.txt', mode='update')
    
    theta = []
    scidata = hdulist[1]
    PSF = hdulist[1].data[1][2]
    PSF = PSF/np.max(PSF)
    THETA = hdulist[2].data
    for i in THETA:
        theta.append(i[0])
    theta = np.array(theta)

    return (theta, PSF)
