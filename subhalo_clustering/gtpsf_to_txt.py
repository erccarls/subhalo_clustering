import sys
sys.path.append('/home/carlson/pyfits/lib')
import pyfits
import numpy as np

# load emissivity data
hdulist = pyfits.open('./psfFront.txt', mode='update')
hdulist2 = pyfits.open('./psfBack.txt', mode='update')
scidata = hdulist[1]
hdulist.info()

# Read FITS header to get dimensions of emissivity map (pixels)

PSF1 = hdulist[1].data[1][2]
PSF2 = hdulist2[1].data[1][2]
THETA1 = hdulist[2].data
THETA2 = hdulist2[2].data

PSF1 = PSF1/np.sum(PSF1)
PSF2 = PSF2/np.sum(PSF2)

cumulative1 = np.cumsum(PSF1/np.sum(PSF1))
cumulative2 = np.cumsum(PSF2/np.sum(PSF2))


Front_95= THETA1[np.argmin(np.abs(cumulative1-.95))]
Back_95= THETA2[np.argmin(np.abs(cumulative2-.95))]

print Front_95
print Back_95

import matplotlib.pyplot as plt
plt.ylim((0,1))
plt.plot(THETA1,PSF1/np.max(PSF1), label="Front", c='b')
plt.plot(THETA1,cumulative1, label="Front Cumulative", c='b', ls = '--')
plt.plot(THETA2,PSF2/np.max(PSF1),label="Back",c = 'r')
plt.plot(THETA2,cumulative2,label="Back Cumulative",c = 'r', ls = '--')

plt.axvline(Front_95, 0, 0.95,color = 'k',ls='-.')
plt.axvline(Back_95, ymin = 0, ymax = 0.95,color = 'k',ls='-.')
plt.axhline(0.95,color = 'k',ls='-.')

plt.ylabel('Containment Fraction')

plt.xlim((0,.5))
plt.xlabel(r'$\theta[^\circ]$')
plt.legend(loc=7)
plt.show()


#iRange = int(dimen[0])
#jRange = int(dimen[1])
#kRange = int(dimen[2])
#lRange = int(dimen[3])
#print "Ranges are", iRange, jRange, kRange, lRange