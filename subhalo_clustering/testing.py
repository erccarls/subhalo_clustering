import pickle
import MC,MCSTATS,DBSCAN
import numpy as np
import sys, math
from matplotlib.backends.backend_pdf import PdfPages
from functools import partial

# Global parameters
angularSize = 10.0 # box surrounding galactic center
outputSize = 300
numProcs = 10


def testing():
    profile = ('PULSAR',)    
    fileOut = 'PULSARRateMap.pickle'
    #mcSims = MC.RUN_PULSAR(5, fileOut, numPhotons=500,numPulsars=1,angularSize=angularSize2, outputSize=outputSize, mcList='test.pickle',HESS=True)
    #MCSTATS.Plot_MC_Positions(mcSims[0], 'test', angularSize = angularSize2)
    
    #pickle.dump(mcSims, open('test.sim','wb'))
    times1, times2 = [],[]
    nphotons = np.logspace(1, 4, 10)
    nphotons = [48,]
    import time
    
    for i in nphotons:
        #profile = ('PULSAR',)
        #map = MC.Gen_Annihilation_Map(angularSize, 300, profile,fileOut)
        mcSims = MC.RUN_PULSAR(10, fileOut, numPhotons=int(i),numPulsars=3,angularSize=10, outputSize=outputSize, mcList='test.pickle')
        #mcSims = MC.RUN(1, fileOut, numPhotons=int(i),angularSize=angularSize2, outputSize=outputSize, mcList='test.pickle',HESS=False)
        
        
        mcSims = pickle.load(open('test.pickle','r'))
        dbscanResults = MCSTATS.DBSCAN_Compute_Clusters(mcSims, eps=.15, min_samples=3, nCorePoints=3, numAnalyze=0)
        
        print MCSTATS.Cluster_Sigs_BG(dbscanResults)
        


import cProfile
import pycallgraph

pycallgraph.start_trace()
cProfile.run('testing()','profile')
pycallgraph.make_dot_graph('test.png')
