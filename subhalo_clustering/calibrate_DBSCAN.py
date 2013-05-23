#===============================================================================
# calibrate_DBSCAN is intended to find optimal DBSCAN parameters for a given
# simulation. It basically the same as in the original clustering paper
# except that we really just don't want background contamination.
# 
# Author: Eric Carlson erccarls@ucsc.edu
# Date Created: 5/20/2013
#===============================================================================

import numpy as np
from matplotlib import pyplot as plt
import MC2, cPickle as pickle
import MCSTATS

# In our benchmark sims we have 115.78^2 square degrees and a background of 2.4e5 photons
# This gives a mean density of rho=17.9 photons per square degree

# First we would like to determine what our detection threshold is, i.e., what is the minimum
# N_min before we begin picking up background clusters. 
# Setup a simulation with just an isotropic background at the density listed above.
# Naively we would expect N_min = rho* pi*eps^2 + n sqrt(rho*pi*eps^2) where n is the rejection significance

# Let us use eps=.5 deg and take an area A=1000*pi*eps^2 -> 28. deg^2

numSims = 1
angularSize = 28.
numBG = 17.9*angularSize**2.
meanBG = numBG/angularSize**2.
eps = .5
deltaTheta = .5


#===============================================================================
# Calibrate N_min
#===============================================================================
#def RunSim():
#    mc = MC2.MC(numPhotons = 0,AngularSize=angularSize) # For now numPhotons is not used, but rather the flux is just the number of photons
#    mc.AddIsotropicSource(flux=numBG)
#    return mc.RunSingle(theta = deltaTheta)
## Run the MC
#mcSims = [RunSim() for i in range(numSims)]
## Scan Parameters
#nMinList = range(12,30,1)
#n_clus = []
#
## Determine number of detected clusters as a function of nMin.  The expected value at 1,2,3 sigma is 18, 22, 25
#for nMin in nMinList:
#    dbscanResults = MCSTATS.DBSCAN_Compute_Clusters(mcSims, eps=eps, min_samples=nMin)    
#    numClusters = len(dbscanResults[0][0])
#    n_clus.append(numClusters)
#    print 'nMin', nMin, 'numClusters/expected',  numClusters/1000.
#    
#from scipy import stats
## Fit Gaussian
#X = np.array(np.linspace(10,30,40))
#NORM = np.max(n_clus)/np.max(stats.norm.pdf(X, loc=14, scale=3.745))
#dist = NORM*stats.norm.pdf(X, loc=14 + 3.745, scale=3.745)
## Plot
#plt.plot(X, dist)
#plt.plot(nMinList ,n_clus)
#plt.xlabel(r'$N_{min}$')
#plt.ylabel(r'$N_{clusters}$')
#plt.show()


# Now that we have established a background rejection threshold (~mean expected BG + 3 sigma), we can tune the epsilon parameter
# We set up another simulation with an isotropic background (same density) and several clusters scattered at 
# random locations, each with the same number of photons and we will scan the number of photons.

numSims = 4
numClusters = 20


#for j in [50,]:
#    def RunSim2():
#        mc = MC2.MC(numPhotons = 0,AngularSize=angularSize) 
#        mc.AddIsotropicSource(flux=numBG)
#        # Sprinkle 100 point sources around
#        locX, locY = (np.random.rand(2,numClusters)-.5)*angularSize 
#        mc.PointSourceLocX    = locX # random locations
#        mc.PointSourceLocY    = locY # random locations
#        mc.PointSourceFlux    = j*np.ones(numClusters)    # 10 sources with j photons per source
#        return mc.RunSingle(theta = deltaTheta)
#    
#    mcSims = [RunSim2() for i in range(numSims)] # Run simulation
#    
#    ## DEBUGGING
#    #plt.scatter(mcSims[0][0],mcSims[0][1],s=1)
#    #plt.show()
##    dbscanResults = MCSTATS.DBSCAN_Compute_Clusters(mcSims, eps=1.5, min_samples=25,plot=True)
##    sigs = MCSTATS.Cluster_Sigs_BG(dbscanResults,BGPhotons= numBG,angularSize = angularSize,numProcs = 1)
##    print len(sigs[0])
##    plt.hist(sigs[0])
##    plt.show()
#    ## END DEBUGGING
#    
#    # Initialize output lists
#    NMIN,EPS = range(10,60,2),np.linspace(.2,.9,15)
#    SIG,NUM = np.zeros((len(NMIN), len(EPS))),np.zeros((len(NMIN), len(EPS)))
#    count =0
#    # Scan over EPS vs NMIN
#    for eps in range(len(EPS)):
#        for nmin in range(len(NMIN)):
#             # Run DBSCAN
#             dbscanResults = MCSTATS.DBSCAN_Compute_Clusters(mcSims, eps=EPS[eps], min_samples=NMIN[nmin])
#             
#             # Input the average significance and number of clusters found
#             for i in range(numSims):
#                 # Compute Significances of Clusters
#                 sigs = MCSTATS.Cluster_Sigs_BG(dbscanResults,BGPhotons=numBG, angularSize = angularSize,numProcs = 1)[i]
#                 # Get a list of members per cluster
#                 members = [len(cluster) for cluster in dbscanResults[i][0]]
#                 if len(members)==0: sig = 0 
#                 else:
#                     sig = np.average(sigs, weights = members)
#                 SIG[nmin][eps] += sig/float(numSims)
#                 print (1+count)/float(numSims*len(EPS)*len(NMIN)), EPS[eps], NMIN[nmin], len(sigs), sig
#                 count+=1
#                 NUM[nmin][eps] += len(dbscanResults[i][0])/float(numSims)
#                          
#    pickle.dump((NMIN,EPS,SIG,NUM), open('eps_vs_nmin'+str(j) ,'wb'))
#    print 'File Saved to eps_vs_nmin'+str(j)
    

#===============================================================================
# Plot
#===============================================================================
ppc = [50,10,20,30] # photons per cluster
for j in range(len(ppc)):
    plt.subplot(2,2,j)
    
    f = open('eps_vs_nmin'+str(ppc[j]),'rb')
    NMIN,EPS,SIG,NUM = pickle.load(f)
    
    # Contourf 
    CS = plt.contourf(EPS,NMIN,SIG,levels = np.linspace(0,10,26))
    # imshow
#    im = plt.imshow(SIG,aspect='auto',vmin=0,vmax=10.,origin='lower', extent = (EPS[0],EPS[-1],NMIN[0],NMIN[-1]))
#    im.set_interpolation('nearest')
#    cbar = plt.colorbar(im)
    
    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = plt.colorbar(CS)
    #cbar.ax.set_ylabel('Global Significance')
    
    CS = plt.contour(EPS,NMIN,NUM,levels = [1,3,10,15,18,19,20,21,22,23,30,30,40,50],colors='k')
    plt.clabel(CS,colors='k')
    # eps-vs-nmin
    EPS1  = np.linspace(0,1, 50)
    NMIN0 = (meanBG*np.pi*np.square(EPS1)) + 0.*(np.sqrt(meanBG*np.pi * np.square(EPS1)))
    NMIN4 = (meanBG*np.pi*np.square(EPS1)) + 3.*(np.sqrt(meanBG*np.pi * np.square(EPS1)))
    plt.plot(EPS1,NMIN0,label=r'MeanBG', linewidth=4)
    plt.plot(EPS1,NMIN4,label=r'MeanBG + $3\sigma$',c='r',linewidth=4)
    plt.text(.7, 15, 'ppc=' + str(ppc[j]))
    
    plt.xlabel(r'$\epsilon$')
    plt.ylabel(r'$N_{min}$')
    plt.xlim((.2,.9))
    plt.ylim((10,60))
    

plt.legend(loc=4,prop={'size':10})
plt.show()


#===============================================================================
# Find Threshold of cluster reconstruction
#===============================================================================
#numSims= 500
#angularSize=7.5
#numBG = 17.8*(angularSize**2.)
#
#
#mc = MC2.MC(numPhotons = 0,AngularSize=angularSize) 
#mc.AddIsotropicSource(flux=numBG)
## Sprinkle 100 point sources around
#def RunSims3(numPhotons):
#    locX, locY = (np.random.rand(2,1)-.5)*angularSize 
#    mc.PointSourceLocX    = [0.,] # random locations
#    mc.PointSourceLocY    = [0.,] # random locations
#    mc.PointSourceFlux    = numPhotons*np.ones(1)    # 10 sources with j photons per source
#    return mc.RunSingle(theta = deltaTheta)
#
#for nmin in [18,22,26]:
#    COUNT,SIGS,SIGS2,NUM_CLUS = [],[],[],[]
#    x = range(0,100,3)
#    for numPhotons in x:
#        mcSims = [RunSims3(numPhotons) for i in range(numSims)]
#        dbscanResults = MCSTATS.DBSCAN_Compute_Clusters(mcSims, eps=.5, min_samples=nmin,indexing=False)
#        
#        count = 0
#        sigs  = 0
#        sigs2 = 0
#        num_clus = 0
#        for sim in dbscanResults:
#            centroidDists = [] 
#            
#            for cluster in sim[0]:
#                X,Y = np.transpose(cluster[0])
#                centroidDists.append( np.sqrt(np.average(X)**2+np.average(Y)**2) )  
#            idx = np.where(np.array(centroidDists)<.75)[0]
#            idx2 = np.where(np.array(centroidDists)>1.)[0]
#            if (len(idx)>0):
#                count += 1
#                #print centroidDists[idx]
#                s = MCSTATS.Cluster_Sigs_BG([sim,],BGPhotons=numBG, angularSize = angularSize,numProcs = 1)[0]
#                if len(idx)==1:
#                    sigs+=s[idx]
#                else: 
#                    sigs+=np.array(s)[idx][0]
#                
#            if (len(idx2)>0):
#                num_clus += 1
#                s = MCSTATS.Cluster_Sigs_BG([sim,],BGPhotons=numBG, angularSize = angularSize,numProcs = 1)[0]
#                if len(idx2)==1:
#                    sigs2+=s[idx2]
#                else: 
#                    sigs2+=np.array(s)[idx2][0]
#            
#                
#        NUM_CLUS.append(num_clus/float(numSims)/17.)
#        COUNT.append(count/float(numSims))
#        if count > 0:
#            SIGS.append(sigs/float(count))
#        else: SIGS.append(0)
#        if num_clus > 0:
#            SIGS2.append(sigs2/float(num_clus))
#        else: SIGS2.append(0)
#    
#    plt.subplot(3,1,1)
#    plt.plot(x, NUM_CLUS,label = str((nmin-14)/4))
#    plt.legend(loc=4,prop={'size':10})
#    plt.ylabel('normed count > 1 deg ')
#    plt.subplot(3,1,2)
#    plt.plot(x, COUNT,label = str((nmin-14)/4))
#    plt.legend(loc=4,prop={'size':10})
#    plt.ylabel('f with 1 cluster < .75 deg ')
#    plt.subplot(3,1,3)
#    plt.plot(x, SIGS,label='real'+str((nmin-14)/4))
#    plt.plot(x, SIGS2,ls='--',label='spurious '+str((nmin-14)/4))
#    plt.ylabel('sig of clusters')
#    plt.xlabel('# Photons')
#    plt.legend(loc=4,prop={'size':10})
#
#plt.show()




    









