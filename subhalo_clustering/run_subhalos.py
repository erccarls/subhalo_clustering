#===============================================================================
# run_subhalos.py
# Run the MC2 code according to CW draft specs
# Author: Eric Carlson 
# Updated: 03/24/2013
#===============================================================================


import numpy as np
from matplotlibn import pyplot as plt
import MC2

# Global Declarations
totalPhotons = int(1e5)
alpha = 1.5
f_s   = .25
mu_res = 20.  # num photons in resolvable limit
angularSize = 50. # simulate square of this width in deg.
deltaTheta = 0.5 # PSF in deg



# Figure out number of photons (MC takes care of Poisson sampling)
num_Sig = f_s*totalPhotons
num_BG  = totalPhotons - num_Sig
omega = (np.pi/180.*angularSize)**2. # Solid angle of simulation
num_subhalos = angularSize**2./(np.pi*deltaTheta**2) # number of point sources to add to simulation

# Find the normalization and define the distribution of subhalo flux 
A = num_sig/(mu_res^(2.-alpha)/(2.-alpha)) # integrate distribution according to eqn. 3 (note mu_res=S_res*exposure so is absorbed into normalization)

# subhalo flux distribution function is dN/dS = A S^{-\alpha} for constant A determined above
# Actually sample the distribution
s = []
while (num_Sig > 0):
    X,Y = np.random(2,num_subhalos)
    s = int(A*X[np.where(Y <= X**-alpha)[0]]) 





