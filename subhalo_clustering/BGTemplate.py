#===============================================================================
# BGTemplate.py: This module generates 2 different background templates based on
#  the fermi galprop BG model, or sideband extrapolation. 
# Author: Eric Carlson
# Updated: 11-14-2012
#===============================================================================
import numpy as np
import ParseFermi
import pickle,sys
import scipy.ndimage
import matplotlib.pyplot as plt
sys.path.append('./pyfits/lib')
import pyfits

def Generate_BG_Template(outputSize=300, angularSize = 10,  fileOut = 'BGRateMap.pickle' ):
    """ 
    Generates background template from 10-200 gev data excluding the 120-140 GeV band 
    
    
    """
    template = np.zeros((outputSize,outputSize))
    ppd=float(outputSize)/float(angularSize) # pixels per deg
    
    events110 = ParseFermi.Import_File('photons.txt', energyRange = (120000,140000),lonRange=(-5,5),latRange = (-5,5))
    events130 = ParseFermi.Import_File('photons.txt', energyRange = (100000,120000),lonRange=(-5,5),latRange = (-5,5))
    events150 = ParseFermi.Import_File('photons.txt', energyRange = (140000,200000),lonRange=(-5,5),latRange = (-5,5))
    
    for i in range(10000,200001,20000):
        if i == 130000:
            continue
        events = ParseFermi.Import_File('photons.txt', energyRange = (i-10000,i+10000),lonRange=(-5,5),latRange = (-5,5))
        BG = np.zeros((outputSize,outputSize))    
        for j in events:
            xIDX = int(j[1]*ppd+float(outputSize/2))
            yIDX = int(j[2]*ppd+float(outputSize/2))
            BG[yIDX][xIDX] += 1.0
            
        psfDeg = .2+float(200)/float(i)
        psfOut = psfDeg*ppd
        #print i/1e3, psfDeg, psfOut
            
        template += scipy.ndimage.filters.gaussian_filter(BG, psfOut)
        
    template = template/np.max(template)
    
    # Write to file        
    outFile = open(fileOut, "wb" )
    pickle.dump(template, outFile)
    print 'Rate Map saved to ', fileOut
    
    plt.imshow(scipy.fliplr(template), 'jet',extent=[5,-5,-5,5])

    plt.xlabel(r'$l [^\circ]$')
    plt.ylabel(r'$b [^\circ]$')
    plt.xlim(5,-5)
    plt.ylim(-5,5)
    plt.colorbar()

    x,y = Find_Centroid(template)
    x,y = (x/ppd -angularSize/2.0,) ,(y/ppd -angularSize/2.0,)
    print x,y
    plt.scatter(x,y, s=10, c='r', marker = '+')
    
    X,Y = FormatEvents(events110)
    plt.scatter(X, Y, label = '100-120 GeV', marker = 'o' , c = 'k')
    
    X,Y = FormatEvents(events130)
    plt.scatter(X, Y, label = '120-140 GeV', marker = 'o' , c = 'r')
    
    X,Y = FormatEvents(events150)
    plt.scatter(X, Y, label = '140-200 GeV', marker = 'o' , c = 'g' )
    
    from matplotlib.font_manager import FontProperties
    fontP = FontProperties()
    fontP.set_size('small')
    plt.legend(loc=1, ncol=1, fancybox=True, shadow=False,prop=fontP,borderaxespad=0.,labelspacing = .2)
    
    from matplotlib.backends.backend_pdf import PdfPages
    if fileOut != '':
        pp = PdfPages(fileOut + '_sideband.pdf')
        plt.savefig(pp, format='pdf')
        print "Figures saved to ", str(fileOut)+ '_sideband.pdf\n',
        pp.close()
    
    plt.show()
    return template



def Generate_BG_Template_Galprop(outputSize=300, angularSize = 10,  fileOut = 'BGRateMap.pickle' ):
    
#    hdulist2 = pyfits.open('gll_psc_v08.fit', mode='update')
#    hdulist2.info()
#    scidata2 = hdulist2[1].data
#    for i in scidata2:
#        if -angularSize/2.0<i[3]<angularSize and -angularSize/2.0<i[4]<angularSize:
#            print i
#    return
    
    ppd=float(outputSize)/float(angularSize) # pixels per deg
    
    events110 = ParseFermi.Import_File('photons.txt', energyRange = (120000,140000),lonRange=(-5,5),latRange = (-5,5))
    events130 = ParseFermi.Import_File('photons.txt', energyRange = (100000,120000),lonRange=(-5,5),latRange = (-5,5))
    events150 = ParseFermi.Import_File('photons.txt', energyRange = (140000,200000),lonRange=(-5,5),latRange = (-5,5))
    
    hdulist = pyfits.open('./gal_2yearp7v6_v0.fits', mode='update')
    hdulist2 = pyfits.open('gll_psc_v08.fit', mode='update')
    scidata = hdulist[0].data
    hdulist.info()
    dimen = scidata.shape
    iRange = int(dimen[0])
    jRange = int(dimen[1])
    kRange = int(dimen[2])
    print "Ranges are", iRange, jRange, kRange
    
    BG = np.zeros((outputSize,outputSize))    
    
    lmid = (kRange-1.0)/2.0
    bmid = (jRange-1.0)/2.0
    scale = 360.0/float(kRange) # Degrees per pixel in the input BG map
    
    for i in range(len(BG[0])):
        for j in range(len(BG[0])):
            XIDX = round(  (i-float(outputSize)/2.0)/float(outputSize)*angularSize/scale + lmid)
            YIDX = round(  (j-float(outputSize)/2.0)/float(outputSize)*angularSize/scale + bmid)
            BG[j][i] = scidata[24][YIDX][XIDX]
    
    template = BG    

    template = template/np.max(template)
    
    # Write to file        
    outFile = open(fileOut, "wb" )
    pickle.dump(template, outFile)
    print 'Rate Map saved to ', fileOut

    plt.imshow(scipy.fliplr(template), 'jet',extent=[5,-5,-5,5])

    plt.xlabel(r'$l [^\circ]$')
    plt.ylabel(r'$b [^\circ]$')
    plt.xlim(5,-5)
    plt.ylim(-5,5)
    plt.colorbar()

    x,y = Find_Centroid(template)
    x,y = (x/ppd -angularSize/2.0,) ,(y/ppd -angularSize/2.0,)
    print x,y
    plt.scatter(x,y, s=10, c='r', marker = '+')
    
    X,Y = FormatEvents(events110)
    plt.scatter(X, Y, label = '100-120 GeV', marker = 'o' , c = 'k')
    
    X,Y = FormatEvents(events130)
    plt.scatter(X, Y, label = '120-140 GeV', marker = 'o' , c = 'r')
    
    X,Y = FormatEvents(events150)
    plt.scatter(X, Y, label = '140-200 GeV', marker = 'o' , c = 'g' )
    
    from matplotlib.font_manager import FontProperties
    fontP = FontProperties()
    fontP.set_size('small')
    plt.legend(loc=1, ncol=1, fancybox=True, shadow=False,prop=fontP,borderaxespad=0.,labelspacing = .2)
    
    from matplotlib.backends.backend_pdf import PdfPages
    if fileOut != '':
        pp = PdfPages(fileOut + '_galprop.pdf')
        plt.savefig(pp, format='pdf')
        print "Figures saved to ", str(fileOut)+ '_galprop.pdf\n',
        pp.close()
    
    plt.show()
    return template



def FormatEvents(events):
    X,Y = [],[]
    for j in events:
        X.append(j[1])
        Y.append(j[2])
    return X,Y


def Find_Centroid(map):
    mass = np.sum(map)
    x,y = 0,0
    
    for i in range(len(map)):
        for j in range(len(map)):
            x += i * map[i][j]/mass
            y += j * map[i][j]/mass
    return x,y
    

#Generate_BG_Template()
#Generate_BG_Template_Galprop()
#Generate_BG_Template_Galprop(angularSize = 4.0,fileOut = 'BGRateMap_HESS_2_deg.pickle')




    
    
    