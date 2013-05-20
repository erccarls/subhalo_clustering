#===============================================================================
# run_subhalos.py
# Run the MC2 code according to CW draft specs
# Author: Eric Carlson 
# Updated: 03/24/2013
#===============================================================================


import numpy as np
from matplotlib import pyplot as plt
import MC2, cPickle as pickle
import multiprocessing as mp
import MCSTATS




#===============================================================================
# Constants
#===============================================================================
deg2rad = np.pi/180.


#===============================================================================
# Configuration
#===============================================================================
numSims      = int(1e0) # Number of simulations to run
N_tot = int(3.2e5) # Total number of signal + BG photons
alpha  = 1.5     # N_eff = (2-alpha) f_s N_max
f_s    = .25     # signal fraction
exposure = 5e10  # exposure in (cm^2 s)
S_res  = 4e-10   # Resolvable flux limit in (cm^2 s)^-1   
mu_res = S_res*exposure  # num photons in resolvable limit
angularSizePercent = 0.325 # percent of sky covered by Fermi
deltaTheta  =  0.5 # 68% containment PSF in deg (using gaussian for now)

angularSize = np.sqrt(4.*np.pi*angularSizePercent / deg2rad**2.) # Simulate square with angular size matching exposure area
#angularSize = 20. # simulate square of this width in deg.
Omega = (deg2rad*angularSize)**2. # Solid angle of simulation


#===============================================================================
# Figure out number of photons and unresolved point sources
#===============================================================================
numSig = f_s*N_tot                  # Number of signal photons
numBG  = N_tot - numSig             # Number of background photons
N_max = N_tot/mu_res                # Maximal number of subhalos (cf. eqn 5)
N_eff = int( (2-alpha)*f_s*N_max)   # Effective number of subhalos

print N_tot, numSig
print N_max, N_eff


import scipy.integrate as integrate
print 

A = numSig/ N_eff / integrate.quad(lambda n: n**(-alpha+1), 0, mu_res)[0] # Normalize distribution according to eqn 3

print 'integrated', N_eff * integrate.quad(lambda x: x**(-alpha)*x,0,20)[0]

maxY = A*(.01**-alpha)
print maxY 
#plt.loglog(np.logspace(-1, 2, 50),A*np.power(np.logspace(-1, 2, 50),-alpha))
#plt.show()
 
 #==============================================================================
 # This runs each monte-carlo.  It is fairly slow because there are so many points
 # in a large angular region 
 #==============================================================================
def Run():
        #===============================================================================
        # This portion happens for each new MC
        #===============================================================================
        # subhalo flux distribution function is dN/dS = A S^{-\alpha}
        # normalization set so that dN/dS=1 when S = 0.1 which should provide good enough sampling
        # Now we sample the distribution
        s = []
        while (N_eff != len(s)):
            X,Y = np.random.rand(2,N_eff-len(s)) # generate random numbers [0,1]
            Y = Y*maxY # max number of photons
            X = X*10*mu_res # rescale the possible fluxes up to max number of photons before being resolved 
            s += list(X[np.where(Y <= A*(X**-alpha))[0]]) # append num photons where Monte Carlo was sampled successfully 
        
        plt.loglog(np.logspace(-1,2,50),A*(np.logspace(-1,2,50)**-alpha)) 
        plt.hist(s,histtype='step')
        plt.show()
        
        print len(s)
        # Assign coordinates (in degrees) to each subhalo point source 
        locX, locY = (np.random.rand(2,N_eff)-.5)*angularSize
        
        # Now we have a number of photons, s[i] assigned to each subhalo i with coordinates locX, locY
        # Initialize the monte-carlo, and add subhalos with the appropriate flux.
        mc = MC2.MC(numPhotons = 0,AngularSize=angularSize) # For now numPhotons is not used, but rather the flux is just the number of photons
        mc.PointSourceLocX    = locX
        mc.PointSourceLocY    = locY
        mc.PointSourceFlux    = s    # list with mean num photons for each source
         
        # Add isotropic BG
        mc.AddIsotropicSource(flux=numBG)
        # Run simulation and return list of coordinate pairs of each photon
        return mc.RunSingle(theta = deltaTheta)


# Initialize and run monte-carlo simulations 
mcSims = [Run() for i in range(numSims)]
#print len(mcSims[0][0])
#plt.scatter(mcSims[0][0], mcSims[0][1], s = 1.)
#plt.show()
# Serialize to file
#pickle.dump(mcSims, open('sims.pickle', 'wb'))



#===============================================================================
# DBSCAN Section
#===============================================================================
eps = .5 # DBSCAN eps parameter
meanBG = numBG / angularSize**2.# Mean BG Density per square deg.  
nMin = np.pi*eps**2.*meanBG + 3*np.sqrt(np.pi*eps**2.*meanBG) # 3 sigma BG fluctuation

#mcSims = pickle.load(open('sims.pickle', 'r'))

#dbscanResults = MCSTATS.DBSCAN_Compute_Clusters(mcSims, eps=eps, min_samples=nMin)
# Return a list of cluster significances for each cluster in each simulation assuming a uniform background contribution
#sigs = MCSTATS.Cluster_Sigs_BG(dbscanResults, BGTemplate = 'BGRateMap.pickle',angularSize = angularSize,BG= 0.75,numProcs = 1)
print sigs
# If we want we could concatenate into a giant list of all clusters
s = np.array(sigs).flatten()
print s
plt.hist(s, 50,histtype='step')
plt.show()








