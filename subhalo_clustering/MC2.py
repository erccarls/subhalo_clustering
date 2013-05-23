#===============================================================================
# MC2.py
# Primary Monte Carlo code which generates photon event maps based on DM subhalo
# models in the MW.  
# Author: Eric Carlson
# Updated: 03/24/2013
#===============================================================================



import numpy as np
import time, pickle
import matplotlib.pyplot as plt
import multiprocessing as mp
from multiprocessing import pool
import FermiPSF

def BuildSubhaloProfile(alpha):
    """
    Returns an radial profile of NFW annihilation projection with rs = 1 and 100 elements.  It is normalized to unit integral
    """
    X = np.linspace(0, 2, 200)
    Y2 = np.square(np.linspace(0, 2, 200))
    # Vectorized profile evaluation.  Factor 10x speedup over pure python.
    rhoDMVec = np.vectorize(rho_DM_NFW)        
    R = np.vectorize(lambda x : np.sum(np.square(rhoDMVec(np.sqrt(x*x+Y2),1.,alpha)) ))(X)
    
    return R/np.max(R)
        
#===============================================================================
# DM Profiles
#===============================================================================
def rho_DM_NFW(r,rc,alpha):
    if (r != 0.0): return (rc/r)**alpha * (1+ r/rc)**(-3+alpha)
    else: return 0.0
def rho_DM_EIN(r,rc,alpha):
    return np.exp(-2.0/alpha*((r/rc)**alpha-1.))


class MC():
    def __init__(self, numPhotons, AngularSize):
        self.AngularSize       = AngularSize        # Size in degrees (assumed cartesian)
        
        self.PointSourceLocX   = []
        self.PointSourceLocY   = []
        self.PointSourceFlux   = []
        
        self.DiffuseSources    = []
        self.NFWSubhalos       = []
        self.IsotropicSources  = []
        self.numPhotons = numPhotons
    
#    class __PointSource():
#        def __init__(self,flux,pos):
#            self.flux = flux # Flux
#            self.pos  = pos # Position
            
    class __NFWSubhalo():
        def __init__(self,flux, pos, rs, alpha):
            self.flux    = flux
            self.pos     = pos
            self.rs      = rs # in degrees because this is the relevent dimension.
            self.alpha   = alpha
            self.profile = BuildSubhaloProfile(alpha) # values of flux 
            self.XX      = np.linspace(0, 2., 200)
                                
    class __IsotropicBG():
        def __init__(self, flux):
            self.flux    = flux
            
    def AddNFWSubhalo(self,flux, pos, rs, alpha):
        # Assign random position
        if (pos ==None): pos = (np.random.rand(2)-.5)*self.AngularSize
        self.NFWSubhalos.append(self.__NFWSubhalo(flux,pos,rs,alpha))
    
    def AddIsotropicSource(self,flux):
        self.IsotropicSources.append(self.__IsotropicBG(flux))
        
#    def AddPointSource(self, flux, pos=None):
#        if (pos ==None): pos = (np.random.rand(2)-.5)*self.AngularSize
#        self.PointSources.append(self.__PointSource(flux,pos))    
    
    #===============================================================================
    # Monte Carlo
    #===============================================================================    
    def __RunSubhalos(self,XMASTER,YMASTER):
        # For each subhalo in the Queue
        XMASTER, YMASTER= [],[]
        for subhalo in self.NFWSubhalos:
            # Choose number of photons according to poisson statistics
            numPhotons = np.random.poisson(subhalo.flux)# Need to scale correctly
            r = []
            # Here we sample from the LOS projected NFW^2 profile
            while (numPhotons != 0):
                X,Y = np.random.rand(2,numPhotons)
                new = list(2.*X[np.where((Y<=np.interp(X*2, subhalo.XX, subhalo.profile)))[0]])
                r+=new
                numPhotons-=len(new)
            # Assign random angles and modify X,Y
            theta = np.random.rand(len(r))*2.*np.pi
            X = subhalo.pos[0] + subhalo.rs*np.cos(theta)*np.array(r)
            Y = subhalo.pos[1] + subhalo.rs*np.sin(theta)*np.array(r)
            # Append to master list
            XMASTER+=list(X)
            YMASTER+=list(Y)
        return XMASTER,YMASTER
    
    def __RunPointSources(self,XMASTER,YMASTER):
        #numPhotons = np.random.poisson(self.PointSourceFlux) # pick num photons for each source
        numPhotons = self.PointSourceFlux # pick num photons for each source
        numPhotonsCum = np.append(np.array([0,]), np.cumsum(numPhotons)) # find cumulative sum
        #print "Num signal photons after poisson sample: ", np.sum(numPhotons)
        XMASTER = np.ones(np.sum(numPhotons)) # create empty arrays (faster than appending lists)
        YMASTER = np.ones(np.sum(numPhotons)) # create empty arrays
        # For all sources, generate a list of photons with initial positions given by the corresponding point source position
        # We will then alter these positions by the PSF below.
        for i in range(len(numPhotons)):
            XMASTER[numPhotonsCum[i]:numPhotonsCum[i+1]] = self.PointSourceLocX[i]   
            YMASTER[numPhotonsCum[i]:numPhotonsCum[i+1]] = self.PointSourceLocY[i]
        
        return XMASTER, YMASTER
    
    def __RunIsotropicSources(self,XMASTER,YMASTER):
        for source in self.IsotropicSources:
            # Choose poisson number of photons
            numPhotons = np.random.poisson((source.flux))
            X,Y = (np.random.rand(2,numPhotons)-.5)*self.AngularSize
            XMASTER += list(X)
            YMASTER += list(Y)
        return XMASTER, YMASTER
            
    def __ApplyPSF(self,XMASTER,YMASTER, psffront,psfback):
        """
            X,Y: List of X and Y vectors to shift
            psffront, psfback are the PSF tables (theta, PSF) from the PSF FITS file
        """
        
        PSFTHETA_FRONT,PSFVALUE_FRONT = psffront
        PSFTHETA_BACK,PSFVALUE_BACK = psfback
        # Find 99.9% containment and don't sample beyond this
        radialCutoff = PSFTHETA_BACK[np.where((np.cumsum(PSFVALUE_BACK)/np.sum(PSFVALUE_BACK) >.999))[0][0]]
        
        FrontArea = 0.561 # Fraction of front converting events.  The rest are back converting.
        numFront  = int(len(XMASTER)*FrontArea) # Number of front converting events
        numBack   = len(XMASTER)-numFront  # Number of back converting events 
        
        r = []
        # Sample Front Events
        while (numFront != 0):
            X,Y = np.random.rand(2,numFront)
            X = radialCutoff*X
            new = list(X[np.where((Y<=np.interp(X, PSFTHETA_FRONT, PSFVALUE_FRONT)))[0]])
            r+=new
            numFront-=len(new)
        # Sample Back Events
        while (numBack != 0):
            X,Y = np.random.rand(2,numBack)
            new = list(X[np.where((Y<=np.interp(X, PSFTHETA_BACK, PSFVALUE_BACK)))[0]])
            r+=new
            numBack-=len(new)
        # Mix up front and back converted events.
        np.random.shuffle(r)
        # Choose Random Angles and Shift Photon Positions
        theta = np.random.rand(len(r))*2.*np.pi
        XMASTER += np.cos(theta)*np.array(r)
        YMASTER += np.sin(theta)*np.array(r)
        return list(XMASTER),list(YMASTER)
    
    def __ApplyGaussianPSF(self, XMASTER, YMASTER, r68):
        r = []
        num = len(XMASTER)
        # given r68, calculate what sigma should be for a gaussian PSF to meet this (can solve analytically)
        sigma = r68/1.0674
        
        while (num != 0):
            X,Y = np.random.rand(2,num)
            X = 5.*sigma*X
            # Normalize to max=1 and perform MC sampling
            new = list(X[np.where(Y <=  np.sqrt(2*np.e)/sigma* X*np.exp((-np.square(X)/(sigma*sigma))  ))[0]])
            r+=new
            num-=len(new)
        
        theta = np.random.rand(len(r))*2.*np.pi
        XMASTER += np.cos(theta)*np.array(r)
        YMASTER += np.sin(theta)*np.array(r)
        
        return list(XMASTER),list(YMASTER)
    
    
    
    #=======================================================================
    # Define single variable function for each simulation
    def RunSingle(self,theta = 0.5):
        np.random.seed() # New random seed
        X,Y= [], []      # Initialize Photon List    
        #X,Y = self.__RunSubhalos(X,Y) # Run Subhalo simulations
        
        start = time.time()
        X,Y = self.__RunPointSources(X,Y) # Run Subhalo simulations
        #print "Finished RunPointSources sources in : " , time.time() - start, ' s'
        
        #X,Y =self.__ApplyPSF(X, Y, PSFTableFront, PSFTableBack) # PSF modulation
        start = time.time()
        X, Y = self.__ApplyGaussianPSF(X,Y, theta)
        #print "Finished apply PSF in : " , time.time() - start, ' s'
        
        # Simulate Isotropic Backgrounds (don't bother with PSF for these)
        start = time.time()
        X,Y = self.__RunIsotropicSources(X,Y)
        #print "Finished isotropic sources in : " , time.time() - start, ' s'
        return X,Y
    #=======================================================================
    
    
    
    
    
    def RunAll(self,numSims,numProcs=1, theta= 0.5):
        """
        Runs All Queued Monte-Carlo Simulations.
        Inputs:
            numSims: number of simulations to use.
            theta: Gaussian PSF in degrees
            numProcs: Number of simultaneous threads to use.
        """
        # Initialize the thread pool
        if (numProcs<=0):numProcs += mp.cpu_count()
        p = pool.Pool(numProcs)
        # Run stats
        #print "Running " , numSims, " simulations using ", numProcs, " thread(s)..." 
        #start = time.time()
        # Load PSF
        PSFTableFront = FermiPSF.PSF_130(convType='front') # load PSF front converting
        PSFTableBack = FermiPSF.PSF_130(convType='back')   # load PSF back converting #Currently 130
        
        #=======================================================================
        # Define single variable function for each simulation
        def RunSingle(i):
            np.random.seed() # New random seed
            X,Y= [], []      # Initialize Photon List    
            #X,Y = self.__RunSubhalos(X,Y) # Run Subhalo simulations
            X,Y = self.__RunPointSources(X,Y) # Run Subhalo simulations
            #X,Y =self.__ApplyPSF(X, Y, PSFTableFront, PSFTableBack) # PSF modulation
            X, Y = self.__ApplyGaussianPSF(X,Y, theta)
            # Simulate Isotropic Backgrounds (don't bother with PSF for these)
            X,Y = self.__RunIsotropicSources(X,Y)
            return X,Y
        #=======================================================================
        
        #result = p.map(RunSingle,range(numSims)) # Multithreaded map
        result = map(RunSingle,range(numSims))    # Serial Map 
        # Run Stats 
        #print 'Ran ', numSims, " simulations in ", time.time()-start , " seconds. (", numSims/(time.time()-start) , " sims/sec)"
        
        return result






        
#mc = MC(numPhotons = 0,AngularSize=20.)

#for i in range(0,10):
    #mc.AddNFWSubhalo(flux=20,pos=None,rs=np.random.rand()*10,alpha=1)
    #mc.AddPointSource(flux=20,pos=None)
#mc.AddIsotropicSource(300)
#X,Y = mc.RunAll(numSims = 1, theta = 0.5, numProcs = 1)[0]
#
#import MCSTATS
#start = time.time()
#print len(MCSTATS.DBSCAN_Compute_Clusters([(X,Y),], .15, min_samples=10, indexing = False)[0][0])
#print time.time()-start
#plt.scatter(X, Y, s=.5)
#plt.show()



#import MCSTATS
#numphotons,time1, time2 = np.logspace(1,6,30), [],[]
#for i in numphotons:
#    mc = MC(numPhotons = 0,AngularSize=50.)
#    mc.AddIsotropicSource(i)
#    X,Y = mc.RunAll(numSims = 1, numProcs = 1)[0]
#    start = time.time()
#    if i<30000:
#        len(MCSTATS.DBSCAN_Compute_Clusters([(X,Y),], .15, min_samples=.15**2*50**2/i+2*np.sqrt(i), indexing = False)[0][0])
#        time1.append(time.time()-start)
#    start = time.time()
#    len(MCSTATS.DBSCAN_Compute_Clusters([(X,Y),], .15, min_samples=.15**2*50**2/i+2*np.sqrt(i), indexing = True)[0][0])
#    time2.append(time.time()-start)
#pickle.dump((numphotons, time1,time2), open('times.pickle','wb'))


#
#numphotons, time1,time2 = pickle.load((open("times.pickle")))
#print np.log10(numphotons)[-3:]
#print time2[-3:]
#plt.plot(numphotons[:len(time1)], time1,label = 'Old')                                        
#plt.plot(numphotons, time2, label='New')
#plt.xlabel('N Photons')
#plt.ylabel('CPU Seconds')
#plt.legend()
#plt.yscale('log')
#plt.xscale('log')
#plt.show()