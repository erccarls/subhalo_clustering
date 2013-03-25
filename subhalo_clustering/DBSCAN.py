#===============================================================================
# DBSCAN.py: Methods for running clustering algorithm and computing some cluster
#  statistics such as significance, background counts, and clustering scale.
# Author: Eric Carlson
# Updated: 03-01-2013
#===============================================================================
import numpy as np
from scipy.spatial import distance
from sklearn.cluster import DBSCAN
import dbscan_indexed
from sklearn import metrics
import pickle, math
from scipy import weave
from scipy.weave import converters

def RunDBScan(X,eps,n_samples,nCorePoints =3 ,indexing = None):
    """
    Runs DBScan on Number-observations-length vector X of coordinate pairs using parameters eps and n_samples defining the search radius and the minimum number of events.
    
    Inputs:
        X:         a list of coordinate pairs (n x 2) for each photon
        eps:       DBSCAN epsilon parameter
        n_samples: DBSCAN core point threshold.  Must have AT LEAST n_samples points per eps-neighborhood to be core point
    Optional inputs 
        nCorePoints: Number of core points required to form a cluster
        indexing:    None: Automatically choose based on number of events
                     True: Grid based indexing (Generally faster and *MUCH* more memory efficient) 
                     False: No indexing.  Computes entire distance matrix.  Faster for < 1000 photons
    Returns: 
    (clusterReturn,labels): a tuple of lists of coordinate pairs for points in each cluster found (core points only, not border points, and labels for each point)
    """
    #===========================================================================
    # Compute DBSCAN
    #===========================================================================
    db = dbscan_indexed.DBSCAN(eps, min_samples=n_samples, indexing = indexing).fit(X)
    
    core_samples = db.core_sample_indices_
    labels = db.labels_
    
    # Get the cluster labels for each core point
    coreLabels = [labels[i] for i in core_samples]
        
    # Count clusters with > nCore, core points
    validClusters = [] 
    [validClusters.append(i) if coreLabels.count(i) >= nCorePoints else None for i in set(coreLabels)]
    
    # relabel points that are not in valid clusters as noise.  If you want border points, comment out this line
    labels = [label if label in validClusters else -1 for label in labels]

    # For each cluster build a list of the core sample coordinates
    X = np.asanyarray(X)
    clusterReturn = [X[np.where((labels == cluster))[0]] for cluster in validClusters]
    
    return (clusterReturn, labels)
    


#===============================================================================
# Internal: Integrates the background template
#===============================================================================
def Evaluate_BG_Contribution(x,y,radius, BGTemplate, numBGEvents, flatLevel = 0): 
    """
    # There is an unresolved bug with this code.  DO NOT USE IN CURRENT FORM
    Integrate the background template and return the expected event count.
    
    Inputs:
     -x,y are the centroid of the cluster 
     -radius is the radius of integration in pixels
     -BGTemplate is the pickled background template used.
     -numBGEvents is the total number of expected background events for the entire angular region being considered. (.75 times)
    
    Returns:
    -count: expected number of background events in region.
    """
    #===========================================================================
    # There is an unresolved bug with this code.  DO NOT USE IN CURRENT FORM 
    #===========================================================================
    # Rescale the BG template so that the integral directly gives the event count.
    BGTemplate = np.array(BGTemplate)/float(np.sum(BGTemplate))*(1.0-flatLevel)
    BGTemplate += flatLevel/np.shape(BGTemplate)[0]**2.0  # Add flat Backgronud
    BGTemplate = float(numBGEvents)*BGTemplate
    
    # Specify data types for weave
    size = len(BGTemplate[0])
    radius = int(round(radius))
    x,y = float(x),float(y)
    start = int(-radius-1)
    

    # Integrate annulus
    code = """
        double ret = 0.;
        for (int i= start; i<-start ; i++){
            for (int j= start; j<-start ; j++){
                if ((i*i+j*j <= radius*radius) && ((0<=(i+x)<size) && (0<=(j+y)<size))){
                    ret += BGTemplate((int)(j+y), (int) (i+x));
                }
            }
        }
        return_val = ret;
    """
    return float(weave.inline(code,['radius','BGTemplate','size','x','y','start'], compiler='gcc', type_converters = converters.blitz)) 
            
    
def Compute_Cluster_Significance(X, BGTemplate, totalPhotons,outputSize=300, angularSize = 10.0,flatLevel = 0,BG = .75):
    """
    Takes input list of coordinate pairs (in angular scale) and computes the cluster significance based on a background model.
    
    Inputs:
        -X is a tuple containing a coordinate pair for each point in a cluster.  
        -BGTemplate is the pickled background array used.
        -totalPhotons is the total number of all photons.  This is used to estimate the background count.
        
        -flatLevel: What fraction of background is isotropic? 1 for completely isotropic 
        -BG: Average fraction of total photons that are background.
    returns significance
    """
    # Default to zero significance
    if (len(X)==1):return 0
    
    # Otherwise.......
    x,y = np.transpose(X) # Reformat input
    numBGEvents = totalPhotons*BG # Number of expected background events.  Based on power law extrapolation from 10-300 GeV excluding 120-140 GeV
    ppa = float(outputSize)/float(angularSize) # pixels per degree
    centX,centY = np.mean(x), np.mean(y) # Compute Centroid
    
    # Build list of radii from cluster centroid
    r = [math.sqrt((x[i]-centX)**2 + (y[i]-centY)**2) for i in range(len(x))]
    
    # Sort the list and choose the radius where the cumulative count is >95%
    countIndex = int(math.ceil(0.95*len(r)-1)) 
    clusterRadius = np.sort(r)[countIndex]   # choose the radius at this index 

    # Estimate the background count
    #N_bg = Evaluate_BG_Contribution(centX*ppa+outputSize/2.0,centY*ppa+outputSize/2.0,clusterRadius*ppa,BGTemplate,numBGEvents, flatLevel = flatLevel)
    
    
    # For now just use isotropic density
    BGDensity = BG*totalPhotons / angularSize**2
    N_bg = np.pi * clusterRadius**2. * BGDensity                 
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
    Computes the mean parwise distance matrix standard deviation 
    Inputs: 
     -cluster: Tuple containing a coordinate pair for each cluster point.
    Returns:
        (mean, std) Mean pairwise dist and std
    '''
    d = distance.pdist(cluster)
    return np.mean(d), np.std(d)
                   
                   
    
    

