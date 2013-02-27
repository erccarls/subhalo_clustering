import pickle
import MC,MCSTATS,DBSCAN
import numpy as np
import sys, math
from matplotlib.backends.backend_pdf import PdfPages
from functools import partial

# Global parameters
angularSize = 10.0 # box surrounding galactic center
angularSize2 = 4.0
size = 300  # size of output ratemap
numSims = 100000 # Number of simulations to use.
numSims2 = 2000
numPhotons = 48
numPhotons2 = 2000
outputSize = 300

numProcs = 10
run48 = False   
run2000 = True

from scipy import weave
cos = np.cos
sin = np.sin
arcsin = np.arcsin
sqrt = np.sqrt

def dist(X,Y):
    # Use Haverside Formula http://en.wikipedia.org/wiki/Great-circle_distance
    n = len(X)
    support = "#include <math.h>"
    d = np.zeros(n)
    code = """
    
    for (int i=1;i<n;i++)
    {
        //printf("%F\\n", X[0]*X[i]);
        d[i] = 2.*asin(sqrt(pow(sin( .5*(Y[0]-Y[i])),2.) + cos(Y[0])*cos(Y[i])*pow(sin(.5*(X[0]-X[i])),2.)));
    }
    """
    weave.inline(code,['n','X','Y','d'],support_code = support, libraries = ['m'])
    
    return 
    

def testing():
#    profile = ('PULSAR',)    
#    fileOut = 'PULSARRateMapHESS.pickle'
    #mcSims = MC.RUN_PULSAR(5, fileOut, numPhotons=500,numPulsars=1,angularSize=angularSize2, outputSize=outputSize, mcList='test.pickle',HESS=True)
    #MCSTATS.Plot_MC_Positions(mcSims[0], 'test', angularSize = angularSize2)
    
    #pickle.dump(mcSims, open('test.sim','wb'))
    times1, times2 = [],[]
    nphotons = np.logspace(1, 4, 10)
    nphotons = [10000,]
    import time
    
    for i in nphotons:
        #mcSims = MC.RUN_PULSAR(10, fileOut, numPhotons=int(i),numPulsars=1,angularSize=angularSize2, outputSize=outputSize, mcList='test.pickle',HESS=False)
        mcSims = MC.RUN(5, fileOut, numPhotons=int(i),angularSize=angularSize2, outputSize=outputSize, mcList='test.pickle',HESS=False)
        
        mcSims = pickle.load(open('test.pickle','r'))
        start = time.time()
        dbscanresults = MCSTATS.DBSCAN_Compute_Clusters(mcSims, eps=.1, n_cluster=30, nCore=3, S_cut=2.0, numAnalyze=0,  HESS=True, angularSize=angularSize2, indexing = False)
        elapsed1 = time.time()-start
        start = time.time()
        dbscanresults = MCSTATS.DBSCAN_Compute_Clusters(mcSims, eps=.1, n_cluster=30, nCore=3, S_cut=2.0, numAnalyze=0,  HESS=True, angularSize=angularSize2, indexing = True)
        elapsed2 = time.time()-start
        print 'nPhotons, Times:', i, elapsed1, elapsed2
        times1.append(elapsed1)
        times2.append(elapsed2)
    pickle.dump((nphotons, elapsed1,elapsed2),open('times.out','wb'))

    
#    import time
#    #mcSims = MC.RUN_PULSAR(1, 'PULSARRateMapHESS.pickle', numPhotons=1000,numPulsars=3,angularSize=angularSize2, outputSize=outputSize, mcList='test.pickle',HESS=False)
#    mcSims = pickle.load(open('test.pickle','r'))
#    start = time.time()
#    dbscanresults = MCSTATS.DBSCAN_Compute_Clusters(mcSims, eps=.1, n_cluster=3, nCore=3, S_cut=2.0, numAnalyze=0,  HESS=True, angularSize=angularSize2, indexing = True)
#    print 'time', time.time()-start




from matplotlib import pyplot as plt
nphotons = [1e2,1e3,2e3,4e3,1e4]
old = [.096,.325,25,192,0]
new = [.184,1.637,3.78,11.1,54] 

plt.plot(nphotons,old, label = 'No Indexing')
plt.plot(nphotons,new, label = 'R-tree Indexing')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('time per core (s)')
plt.xlabel('num photons')
plt.legend()
#plt.show()

#    
    #idx = INDEX.Index()

    #X ,Y = np.array(np.zeros(100000)), np.array(np.ones(100000))
    #dist(X,Y)
    
    
    #[2.*arcsin(sqrt(pow(sin( .5*(Y[0]-Y[i])),2.) + cos(Y[0])*cos(Y[i])*pow(sin(.5*(X[0]-X[i])),2.))) for i in range(100000)]
    
        


import cProfile
import pycallgraph
profile = ('PULSAR',)    
fileOut = 'PULSARRateMapHESS.pickle'
#mcSims = MC.RUN_PULSAR(1, fileOut, numPhotons=100000,numPulsars=1,angularSize=angularSize2, outputSize=outputSize, mcList='test.pickle',HESS=True)


pycallgraph.start_trace()
cProfile.run('testing()','profile')
pycallgraph.make_dot_graph('test.png')
