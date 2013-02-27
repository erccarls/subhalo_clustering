#===============================================================================
# DBSCAN.py: Methods for running clustering algorithm and computing some cluster
#  statistics such as significance, background counts, and clustering scale.
# Author: Eric Carlson
# Updated: 11-14-2012
#===============================================================================
import numpy as np
from scipy.spatial import distance
from sklearn.cluster import DBSCAN
import dbscan_indexed
from sklearn import metrics
import pickle, math

def RunDBScan(X,eps,n_samples,nCore =3 ,plot = False,indexing = True):
    """
    Runs DBScan on Number-observations-length vector X of coordinate pairs using parameters eps and n_samples defining the search radius and the minimum number of events.
    
    If plot==true then a count map with color coded clusters is displayed.
    
    Returns: 
    -clusterReturn: a tuple of lists of coordinate pairs for points in each cluster found (core points only, not border points)
    -n_clusters: the number of clusters found
    -noiseCount: the number of points labeled noise (note that only core border points are not counted noise, and not counted as cluster points.)
    """
    ##############################################################################
    # Compute similarities
    #D = distance.squareform(distance.pdist(X))
    #S = 1 - (D / np.max(D))
    #db = DBSCAN(eps=eps, min_samples=n_samples).fit(S)
    #db = dbscan_indexed.DBSCAN(eps=eps, min_samples=n_samples).fit(S)
    ##############################################################################
    # Compute DBSCAN
    db = dbscan_indexed.DBSCAN(eps, min_samples=n_samples, indexing = indexing).fit(X)
    
    core_samples = db.core_sample_indices_
    labels = db.labels_
    
    # Count clusters with > nCore, core points
    coreLabels = []
    validClusters = []
    for i in core_samples:
        coreLabels.append(labels[i])
    for i in range(int(np.max(labels)+1)):
        if coreLabels.count(i)>=nCore:
            validClusters.append(i)
    # Number of clusters with > nCore core points
    n_clusters_ = int(len(validClusters))
    
    # relabel points that are not in valid clusters as noise
    for i in range(len(labels)):
        if (labels[i] not in validClusters):
            labels[i]=-1
    
    ######################################################
    # For each cluster build a list of the core sample coordinates
    clusterReturn = [] # This is a list of coordinate pairs for points in each cluster.
    
    clusterScales = []
    clusterCounts = []
    #X=zip(X[0],X[1])
    # Form a cluster of core points
    for cluster in validClusters:
        coords = []
        for i in core_samples:
            if labels[i] == cluster:
                coords.append(X[i])                
        clusterReturn.append(coords) # Append a list of coordinate pairs for each cluster

    
    ##############################################################################
    # Plot result
    if (plot == True):   
        import pylab as pl
        from itertools import cycle
        pl.close('all')
        pl.figure(1)
        pl.clf()
        
        # Black removed and is used for noise instead.
        colors = cycle('bgrcmybgrcmybgrcmybgrcmy')
        for k, col in zip(set(labels), colors):
            if k == -1:
                # Black used for noise.
                col = 'k'
                markersize = 6
            class_members = [index[0] for index in np.argwhere(labels == k)]
            cluster_core_samples = [index for index in core_samples
                                    if labels[index] == k]
            for index in class_members:
                x = X[index]
                if index in core_samples and k != -1:
                    markersize = 4
                else:
                    markersize = 2
                pl.plot(x[0], x[1], 'o', markerfacecolor=col,
                        markeredgecolor='k', markersize=markersize)
        
        pl.xlim(5,-5)
        pl.ylim(-5,5)
        pl.title('Estimated number of clusters: %d' % n_clusters_)
        pl.show()
    
    noiseCount = list(labels).count(-1)
    # Compute the average clustering scale with scale weight given by number of elements in a cluster
    return clusterReturn, n_clusters_, noiseCount
    


def Evaluate_BG_Contribution(x,y,radius, BGTemplate, numBGEvents, flatLevel = 0): 
    """
    Integrate the background template and return the expected event count.
    
    Inputs:
     -x,y are the centroid of the cluster, and radius is the radius of integration in pixels
     -BGTemplate is the pickled background template used.
     -numBGEvents is the total number of expected background events for the entire angular region being considered. (.75 times
    
    Returns:
    -count: expected number of background events in region.
    """
    # Rescale the BG template so that the integral directly gives the event count.
    BGTemplate = np.array(BGTemplate)/float(np.sum(BGTemplate))*(1.0-flatLevel)
    BGTemplate += flatLevel/np.shape(BGTemplate)[0]**2.0  # Add flat Backgronud
    BGTemplate = float(numBGEvents)*BGTemplate
    
    size = len(BGTemplate[0])
    count =0.0
    radius = int(round(radius))
    for i in range(-radius-1,radius+1):
        for j in range(-radius-1,radius+1):
            if (i**2+j**2 < radius**2) and 0<(i+x)<size and 0<(j+y)<size:
                count += BGTemplate[j+y][i+x]
    return count
            


    
    
def Compute_Cluster_Significance(X, BGTemplate, totalPhotons,outputSize=300, angularSize = 10.0,flatLevel = 0,SNR = .75):
    """
    Takes input list of coordinate pairs (in angular scale) and computes the cluster significance based on a background model.
    
    Inputs:
        -X is a tuple containing a coordinate pair for each point in a cluster.  
        -BGTemplate is the pickled background array used.
        -totalPhotons is the total number of all photons.  This is used to estimate the background count.
    
    returns significance
    """
    x,y = [],[]
    for i in X:
        x.append(i[0])
        y.append(i[1])
    
    numBGEvents = totalPhotons*SNR # Number of expected background events.  Based on power law extrapolation from 10-300 GeV excluding 120-140 GeV
    
    ppa = float(outputSize)/float(angularSize) # pixels per degree
    
    centX,centY = np.mean(x), np.mean(y) # Compute Centroid
    
    # Build list of radii from cluster centroid
    r = []
    for i in range(len(x)):
        r.append(math.sqrt((x[i]-centX)**2 + (y[i]-centY)**2))
    
    countIndex = int(math.ceil(0.95*len(r)-1)) # Sort the list and choose the radius where the cumulative count is >95%
    
    clusterRadius = np.sort(r)[countIndex]   # choose the radius at this index 
    
    ######################################################
    # Estimate the background count
    N_bg = Evaluate_BG_Contribution(centX*ppa+outputSize/2.0,centY*ppa+outputSize/2.0,clusterRadius*ppa,BGTemplate,numBGEvents, flatLevel = flatLevel)
    N_cl = countIndex - N_bg
    
    
    ######################################################
    # Evaluate significance as defined by Li & Ma (1983)
    if N_cl/(N_cl+N_bg)<1e-20 or N_bg/(N_cl+N_bg)<1e-20:
        return 0
    S2 = 2.0*(N_cl*math.log(2.0*N_cl/(N_cl+N_bg))     +      N_bg*math.log(2.0*N_bg/(N_cl+N_bg)))
    if S2>0.0:
        return math.sqrt(S2)   
    else:
        return 0
    
    
def Compute_Cluster_Scale(cluster):
    '''
    Computes the parwise distance matrix and takes the average of the flattened upper triangular matrix. 
    
    Inputs: 
     -cluster: Tuple containing a coordinate pair for each cluster point.
    Returns the cluster scale
    '''
    d = distance.squareform(distance.pdist(cluster))
    distances = []
    
    # Take upper triangular elements and average.  This gives the average distance from on point to another
    for i in range(len(d)):
        for j in range(i+1,len(list(d[i]))):
            distances.append(d[i][j])
    return np.average(distances), np.std(distances)
    

