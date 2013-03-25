#===============================================================================
# run_subhalos.py
# Run the MC2 code according to CW draft specs
# Author: Eric Carlson 
# Updated: 03/24/2013
#===============================================================================


import numpy as np
from matplotlib import pyplot as plt
import MC2, pickle
import multiprocessing as mp
import MCSTATS


# Configureation and Global Declarations
numSims      = int(1e0)
totalPhotons = int(1e3)
alpha  = 1.5  # 
f_s    = .25  # signal fraction
mu_res = 20.  # num photons in resolvable limit
angularSize = 20. # simulate square of this width in deg.
deltaTheta  = 0.5 # PSF in deg





#=========================================================================
# Figure out number of photons (MC takes care of Poisson sampling)
numSig = f_s*totalPhotons
numBG  = totalPhotons - numSig
omega = (np.pi/180.*angularSize)**2. # Solid angle of simulation
numSubhalos = int(angularSize**2./(np.pi*deltaTheta**2)) # number of point sources to add to simulation
# Normalize Subhalo Flux to total number of signal photons  
A = numSig/(mu_res**(2.-alpha) / (2.-alpha)) # integrate distribution according to eqn. 3 (note mu_res=S_res*exposure so is absorbed into normalization)


 #==============================================================================
 # This runs each monte-carlo.  It is fairly slow because there are so many point
 # in a large angular region 
 #==============================================================================
def Run():
        #===============================================================================
        # This portion happens for each new MC
        #===============================================================================
        # subhalo flux distribution function is dN/dS = A S^{-\alpha} for constant A determined above
        # Actually sample the distribution
        s = []
        while (numSubhalos != len(s)):
            X,Y = np.random.rand(2,numSubhalos-len(s)) # gen random numbers
            X = X*mu_res
            s += list(A * X[np.where(Y <= X**-alpha)[0]]) # append list where the monte-carlo sampling was valid
        
        # Now we have a number of photons, s[i] assigned to each subhalo i
        # Initialize the monte-carlo, and add subhalos with the appropriate flux.
        mc = MC2.MC(numPhotons = 0,AngularSize=angularSize) # For now numPhotons is not used, but rather the flux is just the number of photons 
        # Add each point source  NEED TO VECTORIZE THIS STEP
        [mc.AddPointSource(flux=numPhotons,pos=None) for numPhotons in s]
        # Add isotropic BG
        mc.AddIsotropicSource(flux=numBG)
        # Run simulation and return list of coordinate pairs
        return mc.RunSingle(theta = deltaTheta)


# Initialize and run monte-carlo simulations 
mcSims = [Run() for i in range(numSims)]
# Serialize to file
pickle.dump(mcSims, open('sims.pickle', 'wb'))



#===============================================================================
# DBSCAN Section
#===============================================================================
eps = .5 # DBSCAN eps parameter
meanBG = numBG / angularSize**2.# Mean BG Density per square deg.  
nMin = np.pi*eps**2.*meanBG + 3*np.sqrt(np.pi*eps**2.*meanBG) # 3 sigma BG fluctuation

mcSims = pickle.load(open('sims.pickle', 'r'))
dbscanResults = MCSTATS.DBSCAN_Compute_Clusters(mcSims, eps=eps, min_samples=nmin)

# Return a list of cluster significances for each cluster in each simulation assuming a uniform background contribution
sigs = Cluster_Sigs_BG(dbscanResults, BGTemplate = 'BGRateMap.pickle',angularSize = 10.,BG= 0.75,numProcs = 1)

# If we want we could concatenate into a giant list of all clusters
s = []
[s+=list(sig) for sig in sigs]








