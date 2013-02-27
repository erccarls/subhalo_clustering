import matplotlib.pyplot as plt  #@UnresolvedImport
import matplotlib.image as mpimg #@UnresolvedImport
import matplotlib.cm as cm
import matplotlib, scipy
import pickle
import numpy as np
import scipy.cluster as cluster
from matplotlib.backends.backend_pdf import PdfPages
import csv,sys
sys.path.append('/home/carlson/pyfits/lib')
import pyfits

def Import_PSF(CSVFile,energyRange = (128.0,131.0)):
    '''
    Internal function.  Deprecated.
    0. Energy, 
    1. Livetime, 
    2. Radial position from center (degrees), 
    3. PSF value in (sr^-1)
    '''
    values= []
    csvread = csv.reader(open(CSVFile, 'rb'),delimiter = ' ' )
    for i in csvread:
        if (float(i[0])>energyRange[0] and float(i[0])<energyRange[1]):  
            values.append((float(i[0]),float(i[2]),float(i[3])))
    return values

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
    

#===============================================================================
# These are for specific plots and are not general methods
#===============================================================================
#x,y = PSF_130()
#plt.plot(y[:-5500])
#plt.show()


#plt.xlabel(r'r$[^\circ ]$')
#plt.ylabel('PSF')
#
#
#PSF = Import_PSF('psf_energies.txt',energyRange = (100.0,152.0))
#PSF = Import_PSF('psf_energies.txt',energyRange = (10000.0,101000.0))
#fig = plt.figure()
#ax1 = fig.add_subplot(2,1,1)
#ax1.set_ylabel(r'PSF',fontsize=12)
#plt.xlim((0,.25))
#


#
#PSF = Import_PSF('psf_energies.txt',energyRange = (10000.0,101000.0))
#r = []
#psf = []
#psf2 = [] # cumulative
#for i in PSF:
#    if round(i[0])==100000.0:
#        r.append(i[1])
#        psf.append(i[2]*i[1])
#        psf2.append(np.sum(psf))
#psf = np.array(psf)/np.sum(psf)
#psf2 = np.array(psf2)/np.max(psf2)
#
#for i in range(len(psf2)):
#    if psf2[i]>.68:
#        print '68% containment at r=' + str(r[i])
#        break
#
#plt.plot(r, psf)
#plt.xlim(0,.3)
#plt.show()
#
#plt.clf()
#plt.xlim(0,.3)
#plt.plot(r,psf2)
#plt.show()




#r = []
#psf = []
#psf2 = []
#for i in PSF:
#    if round(i[0])==100000.0:
#        r.append(i[1])
#        psf.append(i[2])
#plt.plot(r, psf)
#
#r = []
#psf = []
#psf2 = []
#for i in PSF:
#    if round(i[0])==50003.0:
#        r.append(i[1])
#        psf.append(i[2])
#plt.plot(r, psf)
#print np.max(psf)
#
#r = []
#psf = []
#psf2 = []
#for i in PSF:
#    if round(i[0])==10155.0:
#        r.append(i[1])
#        psf.append(i[2])
#plt.plot(r, psf)
#
#
#plt.legend(('100 GeV','50 GeV','10 GeV'))
#
#r = []
#psf = []
#psf2 = []
#for i in PSF:
#    if round(i[0])==132.0:
#        r.append(i[1])
#        psf.append(i[2])
#plt.plot(r, psf)
#
#r = []
#psf = []
#psf2 = []
#for i in PSF:
#    if round(i[0])==152:
#        r.append(i[1])
#        psf.append(i[2])
#plt.plot(r, psf)
##plt.legend(('102 GeV','126 GeV','129 GeV','132 GeV','150 GeV')
#
#ax2 = fig.add_subplot(2,1,2)
#ax2.set_ylabel(r'PSF*R',fontsize=12)
#ax2.set_xlabel(r'r$[^\circ ]$',fontsize=12)
#plt.xlim((0,.25))
#
#r = []
#psf = []
#psf2 = []
#for i in PSF:
#    if round(i[0])==100000.0:
#        r.append(i[1])
#        psf.append(i[2]*i[1]/105602.35)
#plt.plot(r, psf)
#
#
#r = []
#psf = []
#psf2 = []
#for i in PSF:
#    if round(i[0])==50003.0:
#        r.append(i[1])
#        psf.append(i[2]*i[1]/105602.35)
#plt.plot(r, psf)
#
#r = []
#psf = []
#psf2 = []
#for i in PSF:
#    if round(i[0])==10155.0:
#        r.append(i[1])
#        psf.append(i[2]*i[1]/105602.35)
#plt.plot(r, psf)
#print np.sum(psf)
#r = []
#psf = []
#psf2 = []
#for i in PSF:
#    if round(i[0])==132.0:
#        r.append(i[1])
#        psf.append(i[2]*i[1]/105602.35)
#plt.plot(r, psf)
#print np.sum(psf)
#
#r = []
#psf = []
#psf2 = []
#for i in PSF:
#    if round(i[0])==152:
#        r.append(i[1])
#        psf.append(i[2]*i[1]/105602.35)
#plt.plot(r, psf)
#print np.sum(psf)
#
#fileOut ='Fermi_PSF'
##plt.legend(('102 GeV','126 GeV','129 GeV','132 GeV','150 GeV'))
#pp = PdfPages(fileOut + '.pdf')
#plt.savefig(pp, format='pdf')
#print "Figures saved to ", str(fileOut)+ '.pdf\n',
#pp.close()
#
#plt.show()
