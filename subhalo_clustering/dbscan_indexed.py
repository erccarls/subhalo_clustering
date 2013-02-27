# -*- coding: utf-8 -*-
"""
DBSCAN: Density-Based Spatial Clustering of Applications with Noise

Modified from Sklearn libraries.  Incorporates R-tree spatial indexing
"""


import warnings
import numpy as np
 
from sklearn.base import BaseEstimator, ClusterMixin
from sklearn.metrics import pairwise_distances
from sklearn.utils import check_random_state

from functools import partial
from scipy import weave

def dbscan2(X, eps=0.5, min_samples=5, metric='euclidean'):
    """Perform DBSCAN clustering from vector array or distance matrix.

    Parameters
    ----------
    X: array [n_samples, n_samples]
    eps: float, optional
        The maximum distance between two samples for them to be considered
        as in the same neighborhood.
    min_samples: int, optional
        The number of samples in a neighborhood for a point to be considered
        as a core point.
    metric: string
        Compute distances in euclidean, or spherical coordinate space

    Returns
    -------
    core_samples: array [n_core_samples]
        Indices of core samples.

    labels : array [n_samples]
        Cluster labels for each point.  Noisy samples are given the label -1.

    References
    ----------
    Ester, M., H. P. Kriegel, J. Sander, and X. Xu, “A Density-Based
    Algorithm for Discovering Clusters in Large Spatial Databases with Noise”.
    In: Proceedings of the 2nd International Conference on Knowledge Discovery
    and Data Mining, Portland, OR, AAAI Press, pp. 226–231. 1996
    """
    from rtree import index as INDEX
    # Get data into correct format
    #X = np.asarray(zip(X[0],X[1]))   # convert to tuples of (long,lat)
    X = np.asarray(X)    # convert to numpy array
    n = np.shape(X)[0]   # Number of points
    deg2rad = np.pi/180. # Conversion
    
    ######################################################################   
    # Build r-tree index.  Bulk loading is slower, don't follow advise of
    # R-tree docs.  Either way this step is very fast.
    idx = INDEX.Index()
    # Add each photon to index
    [idx.insert(i,(X[i][0], X[i][1],X[i][0],X[i][1])) for i in range(0,n)]
    
    # Avoid dots in loops for speed!
    intersection = idx.intersection
    ########################################################################
    # Now that index is built we can simply query the index by calling the 
    # following function.  Efficiently returns square around a point. If in
    # spherical coordinates, we must scale the longitude by 1/cos(latitude)
    # Later we will refine this epsilon neighborhood by actually calculating 
    # distances, so don't need to be too careful here.
    def queryEps(pos,eps):
        if (metric == 'euclidean'):
            return list(intersection((pos[0]-eps, pos[1]-eps,pos[0]+eps,pos[1]+eps)))
        elif (metric == 'spherical'):
            # Be careful near poles, just return larger box in this case.
            if(abs(90.-pos[1]) < eps):
                return list(intersection((-180., 90.-2*eps,180.,90)))
            elif (abs(90.+pos[1]) < eps):
                return list(intersection((-180., -90.,180.,-90+2*eps)))
            # rescale the longitude bounds depending on the latitude
            eps2 = eps/abs(np.cos(pos[1])) 
            return list(intersection((pos[0]-eps2, pos[1]-eps,pos[0]+eps2,pos[1]+eps)))
        else: 
            print "Invalid DBSCAN Metric."
            return None
        
    
    ############################################################################
    # This section loops through the rtree neighborhood and checks if points are
    # really in neighborhood.  Note R-tree returns a rectangular window so we will
    # always need to do this without modifying the rtree package.  It also allows
    # for precise computation of the central angle in spherical coordinates.  
    def refine_nhood(nhood,index,eps):
        if (metric == 'euclidean'):
            x1, y1 = float(X[index][0]), float(X[index][1])
            x2 = np.asarray([X[i][0] for i in nhood])
            y2 = np.asarray([X[i][1] for i in nhood])
            n2 = len(x2)
            support = "#include <math.h>"
            code = """
            py::list ret;
            for (int i=0;i<n2;i++){
                double d = sqrt( (y1-y2[i])*(y1-y2[i]) + (x1-x2[i])*(x1-x2[i]) );
                if (d <= eps) ret.append(nhood[i]);
            }
            return_val = ret;
            """
            return np.array(weave.inline(code,['n2','x2','y2','nhood','x1','y1','eps'],support_code = support, libraries = ['m']))    
            
            
            #def neighbor(p): return (np.sqrt((X[index][0]-X[p][0])*(X[index][0]-X[p][0]) + (X[index][1]-X[p][1])*(X[index][1]-X[p][1]))<=eps) 
            #return np.array(filter(neighbor,nhood))
            
        if (metric == 'spherical'):
            # Use Haverside Formula http://en.wikipedia.org/wiki/Great-circle_distance
            x1, y1 = float(X[index][0]), float(X[index][1])
            x2 = np.asarray([X[i][0] for i in nhood])
            y2 = np.asarray([X[i][1] for i in nhood])
            n2 = len(x2)
            support = "#include <math.h>"
            epsRad = eps*deg2rad
            code = """
            py::list ret;
            for (int i=0;i<n2;i++){
                double d = 2.*asin(sqrt(pow(sin( deg2rad*.5*(y1-y2[i])),2.) + cos(deg2rad*y1)*cos(deg2rad*y2[i])*pow(sin(deg2rad*.5*(x1-x2[i])),2.)));
                if (d <= epsRad) ret.append(nhood[i]);
            }
            return_val = ret;
            """
            return np.array(weave.inline(code,['n2','x2','y2','nhood','x1','y1','epsRad','deg2rad'],support_code = support, libraries = ['m']))    
    
    
    
    # Get rough neighborhoods
    neighborhoods = [queryEps(X[i], eps) for i in range(0,n)]  # get a subspace using spatial index
    
    import time
    start = time.time()
    
    
    # Refine neighborhoods
    neighborhoods = [refine_nhood(neighborhoods[index],index,eps) for index in range(0,n)] # compute dist matrix in subspace
    print 'elapsed' ,time.time()-start
    
    
    
    #======================================================
    # From here the algorithm is essentially the same
    #======================================================
    # Initially, all samples are noise.
    labels = -np.ones(n)
    # A list of all core samples found.
    core_samples = []
    # label_num is the label given to the new cluster
    label_num = 0

    # Look at all samples and determine if they are core.
    # If they are then build a new cluster from them.
    for index in range(0,n):
        if labels[index] != -1 or len(neighborhoods[index]) < min_samples:
            # This point is already classified, or not enough for a core point.
            continue
        core_samples.append(index)

        labels[index] = label_num
        # candidates for new core samples in the cluster.
        candidates = [index]
        while len(candidates) > 0:
            new_candidates = []
            # A candidate is a core point in the current cluster that has
            # not yet been used to expand the current cluster.
            for c in candidates:
                noise = np.where(labels[neighborhoods[c]] == -1)[0]
                noise = neighborhoods[c][noise]
                labels[noise] = label_num
                for neighbor in noise:
                    # check if its a core point as well
                    if len(neighborhoods[neighbor]) >= min_samples:
                        # is new core point
                        new_candidates.append(neighbor)
                        core_samples.append(neighbor)
            # Update candidates for next round of cluster expansion.
            candidates = new_candidates
        # Current cluster finished.
        # Next core point found will start a new cluster.
        label_num += 1
    return core_samples, labels



def dbscan3(X, eps=0.5, min_samples=5, metric='euclidean', indexing =True):
    """Perform DBSCAN clustering from vector array or distance matrix.

    Parameters
    ----------
    X: array [n_samples, n_samples]
    eps: float, optional
        The maximum distance between two samples for them to be considered
        as in the same neighborhood.
    min_samples: int, optional
        The number of samples in a neighborhood for a point to be considered
        as a core point.
    metric: string
        Compute distances in euclidean, or spherical coordinate space

    Returns
    -------
    core_samples: array [n_core_samples]
        Indices of core samples.

    labels : array [n_samples]
        Cluster labels for each point.  Noisy samples are given the label -1.

    References
    ----------
    Ester, M., H. P. Kriegel, J. Sander, and X. Xu, “A Density-Based
    Algorithm for Discovering Clusters in Large Spatial Databases with Noise”.
    In: Proceedings of the 2nd International Conference on Knowledge Discovery
    and Data Mining, Portland, OR, AAAI Press, pp. 226–231. 1996
    """
    
    
    # Get data into correct format
    #X = np.asarray(zip(X[0],X[1]))   # convert to tuples of (long,lat)
    X = np.asarray(X,dtype = np.float32)    # convert to numpy array
    n = np.shape(X)[0]   # Number of points
    deg2rad = np.pi/180. # Conversion
    where = np.where
    sum=np.sum
    square = np.square
    
    #======================================================================
    # Grid based indexing will simply assign each point a set of grid
    # coordinates.  Grid resolution is boxes of width 2*epsilon thus for an
    # epsilon query, we will generally need to query the nearest four boxes.
    # numpy is used to quickly compute grid positions.
    #GRID = np.round(X/(eps))
    # The direction of the four boxes we need to query is determined by the 
    # rounded off part. Positive difference indicates query in the larger
    # index direction
    #GRID_MOD = np.sign(X/eps - GRID)
    
    
    
    ############################################################################
    # This section loops through the rtree neighborhood and checks if points are
    # really in neighborhood.  Note R-tree returns a rectangular window so we will
    # always need to do this without modifying the rtree package.  It also allows
    # for precise computation of the central angle in spherical coordinates.  
    def refine_nhood(nhood,index,eps):
        if (metric == 'euclidean'):
            x1, y1 = float(X[index][0]), float(X[index][1])
            x2 = np.asarray([X[i][0] for i in nhood])
            y2 = np.asarray([X[i][1] for i in nhood])
            n2 = len(x2)
            support = "#include <math.h>"
            code = """
            py::list ret;
            for (int i=0;i<n2;i++){
                double d = sqrt( (y1-y2[i])*(y1-y2[i]) + (x1-x2[i])*(x1-x2[i]) );
                if (d <= eps) ret.append(nhood[i]);
            }
            return_val = ret;
            """
            return np.array(weave.inline(code,['n2','x2','y2','nhood','x1','y1','eps'],support_code = support, libraries = ['m']))    
            
            
            #def neighbor(p): return (np.sqrt((X[index][0]-X[p][0])*(X[index][0]-X[p][0]) + (X[index][1]-X[p][1])*(X[index][1]-X[p][1]))<=eps) 
            #return np.array(filter(neighbor,nhood))
            
        if (metric == 'spherical'):
            # Use Haverside Formula http://en.wikipedia.org/wiki/Great-circle_distance
            x1, y1 = float(X[index][0]), float(X[index][1])
            x2 = np.asarray([X[i][0] for i in nhood])
            y2 = np.asarray([X[i][1] for i in nhood])
            n2 = len(x2)
            support = "#include <math.h>"
            epsRad = eps*deg2rad
            code = """
            py::list ret;
            for (int i=0;i<n2;i++){
                double d = 2.*asin(sqrt(pow(sin( deg2rad*.5*(y1-y2[i])),2.) + cos(deg2rad*y1)*cos(deg2rad*y2[i])*pow(sin(deg2rad*.5*(x1-x2[i])),2.)));
                if (d <= epsRad) ret.append(nhood[i]);
            }
            return_val = ret;
            """
            return np.array(weave.inline(code,['n2','x2','y2','nhood','x1','y1','epsRad','deg2rad'],support_code = support, libraries = ['m']))    
    
    # This line looks for all points with indices in one of the four squares bordering the original point
    # i.e. it builds the initial neighborhoods.  These need to be refined later
    
    #GRID_MOD = GRID+GRID_MOD
    GTX,GTY = np.transpose(X)
    
    #print GRID[0][0]-1
    #print GTX
    
    #======================================================
    # From here the algorithm is essentially the same
    #======================================================
    # Initially, all samples are noise.
    labels = -np.ones(n)

    def epsQuery(i):
        #xcut = np.where( ((GRID[i][0]-2 <= GTX ) & (GTX <= GRID[i][0]+2)))[0]
        #return xcut[np.where(  ((GRID[i][1]-2 <= GTY[xcut]) & (GTY[xcut] < GRID[i][1]+2) ))[0]]
        xcut = where( ((X[i][0]-eps <= GTX ) & (GTX <= X[i][0]+eps)))[0]
        return xcut[where(  ((X[i][1]-eps <= GTY[xcut]) & (GTY[xcut] < X[i][1]+eps) ))[0]]
    neighborhoods = [epsQuery(i) for i in range(0,n)]
    #neighborhoods = np.asarray([where( np.all(np.less_equal(np.abs(GRID_MOD-GRID[i]),[2,2]),axis=1 ))[0] for i in range(0,n)])
    
    ndist = 0
    
    
    #def refine(i):
        #if len(neighborhoods[i])>=min_samples:
            
        #else: return 
    neighborhoods = [(neighborhoods[i][where( sum(square( X[neighborhoods[i]] - X[i]),axis=1) <= eps*eps)[0]]) if len(neighborhoods[i])>=min_samples else neighborhoods[i] for i in range(0,n)]
    # Refine this neighborhood
    #neighborhoods = np.asarray([  neighborhoods[i][where( sum(square( X[neighborhoods[i]] - X[i]),axis=1) <= eps*eps)[0]] for i in range(0,n)  ])    

    # Refine neighborhoods
    #neighborhoods = [refine_nhood(neighborhoods[index],index,eps) for index in range(0,n)] # compute dist matrix in subspace
    
    #======================================================
    # From here the algorithm is essentially the same
    #======================================================
    # A list of all core samples found.
    core_samples = []
    # label_num is the label given to the new cluster
    label_num = 0

    # Look at all samples and determine if they are core.
    # If they are then build a new cluster from them.
    for index in range(0,n):
        if labels[index] != -1 or len(neighborhoods[index]) < min_samples:
            # This point is already classified, or not enough for a core point.
            continue
        core_samples.append(index)

        labels[index] = label_num
        # candidates for new core samples in the cluster.
        candidates = [index]
        while len(candidates) > 0:
            new_candidates = []
            # A candidate is a core point in the current cluster that has
            # not yet been used to expand the current cluster.
            for c in candidates:
                noise = np.where(labels[neighborhoods[c]] == -1)[0]
                noise = neighborhoods[c][noise]
                labels[noise] = label_num
                for neighbor in noise:
                    # check if its a core point as well
                    if len(neighborhoods[neighbor]) >= min_samples:
                        # is new core point
                        new_candidates.append(neighbor)
                        core_samples.append(neighbor)
            # Update candidates for next round of cluster expansion.
            candidates = new_candidates
        # Current cluster finished.
        # Next core point found will start a new cluster.
        label_num += 1
    #print "Core Samples", len(core_samples), " Distance Comps: ", ndist
    return core_samples, labels


def dbscan(X, eps=0.5, min_samples=5, metric='euclidean', indexing = False):
    """Perform DBSCAN clustering from vector array or distance matrix.

    Parameters
    ----------
    X: array [n_samples, n_samples] or [n_samples, n_features]
        Array of distances between samples, or a feature array.
        The array is treated as a feature array unless the metric is given as
        'precomputed'.
    eps: float, optional
        The maximum distance between two samples for them to be considered
        as in the same neighborhood.
    min_samples: int, optional
        The number of samples in a neighborhood for a point to be considered
        as a core point.
    metric: string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string or callable, it must be one of
        the options allowed by metrics.pairwise.calculate_distance for its
        metric parameter.
        If metric is "precomputed", X is assumed to be a distance matrix and
        must be square.
    numCPU: number of cpus.  If < 0. numSystemCPU+1+numCPU are used   
    Returns
    -------
    core_samples: array [n_core_samples]
        Indices of core samples.

    labels : array [n_samples]
        Cluster labels for each point.  Noisy samples are given the label -1.

    Notes
    -----
    See examples/plot_dbscan.py for an example.

    References
    ----------
    Ester, M., H. P. Kriegel, J. Sander, and X. Xu, “A Density-Based
    Algorithm for Discovering Clusters in Large Spatial Databases with Noise”.
    In: Proceedings of the 2nd International Conference on Knowledge Discovery
    and Data Mining, Portland, OR, AAAI Press, pp. 226–231. 1996
    """
    import time
    start = time.time()
    X = np.asarray(X)
    n = X.shape[0]
    # If index order not given, create random order.
    D = pairwise_distances(X, metric=metric)
    
    #from scipy.spatial import distance
    #deg2rad = np.pi/180.
    #def haverside(u,v):
    #    return 2.*np.arccos(np.sqrt(  np.square(np.sin(deg2rad*.5*(u[1]-v[1])) ) + np.cos(deg2rad*u[1])*np.cos(deg2rad*v[1]) * np.square(np.sin(deg2rad*.5 *(u[0]-v[1]))) ))    
    #D = distance.squareform(distance.pdist(X, haverside))
     
    # Calculate neighborhood for all samples. This leaves the original point
    # in, which needs to be considered later (i.e. point i is the
    # neighborhood of point i. While True, its useless information)
    neighborhoods = [np.where(x <= eps)[0] for x in D]
    
    # Initially, all samples are noise.
    labels = -np.ones(n)
    # A list of all core samples found.
    core_samples = []
    # label_num is the label given to the new cluster
    label_num = 0
    # Look at all samples and determine if they are core.
    # If they are then build a new cluster from them.
    for index in range(0,n):
        if labels[index] != -1 or len(neighborhoods[index]) < min_samples:
            # This point is already classified, or not enough for a core point.
            continue
        core_samples.append(index)
        labels[index] = label_num
        # candidates for new core samples in the cluster.
        candidates = [index]
        while len(candidates) > 0:
            new_candidates = []
            # A candidate is a core point in the current cluster that has
            # not yet been used to expand the current cluster.
            for c in candidates:
                noise = np.where(labels[neighborhoods[c]] == -1)[0]
                #print noise, neighborhoods[c][noise]
                noise = neighborhoods[c][noise]
                labels[noise] = label_num
                for neighbor in noise:
                    # check if its a core point as well
                    if len(neighborhoods[neighbor]) >= min_samples:
                        # is new core point
                        new_candidates.append(neighbor)
                        core_samples.append(neighbor)
            # Update candidates for next round of cluster expansion.
            candidates = new_candidates
        # Current cluster finished.
        # Next core point found will start a new cluster.
        label_num += 1
    #print "Core Samples", len(core_samples)," Distance Comps: ", n**2
    return core_samples, labels


class DBSCAN(BaseEstimator, ClusterMixin):
    """Perform DBSCAN clustering from vector array or distance matrix.

    DBSCAN - Density-Based Spatial Clustering of Applications with Noise.
    Finds core samples of high density and expands clusters from them.
    Good for data which contains clusters of similar density.

    Parameters
    ----------
    eps : float, optional
        The maximum distance between two samples for them to be considered
        as in the same neighborhood.
    min_samples : int, optional
        The number of samples in a neighborhood for a point to be considered
        as a core point.
    metric : string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string or callable, it must be one of
        the options allowed by metrics.pairwise.calculate_distance for its
        metric parameter.
        If metric is "precomputed", X is assumed to be a distance matrix and
        must be square.
    random_state : numpy.RandomState, optional
        The generator used to initialize the centers. Defaults to numpy.random.

    Attributes
    ----------
    `core_sample_indices_` : array, shape = [n_core_samples]
        Indices of core samples.

    `components_` : array, shape = [n_core_samples, n_features]
        Copy of each core sample found by training.

    `labels_` : array, shape = [n_samples]
        Cluster labels for each point in the dataset given to fit().
        Noisy samples are given the label -1.

    Notes
    -----
    See examples/plot_dbscan.py for an example.

    References
    ----------
    Ester, M., H. P. Kriegel, J. Sander, and X. Xu, “A Density-Based
    Algorithm for Discovering Clusters in Large Spatial Databases with Noise”.
    In: Proceedings of the 2nd International Conference on Knowledge Discovery
    and Data Mining, Portland, OR, AAAI Press, pp. 226–231. 1996
    """

    def __init__(self, eps=0.5, min_samples=5, metric='euclidean', indexing=None):
        self.eps = eps
        self.min_samples = min_samples
        self.metric = metric
        self.indexing = indexing
        
    def fit(self, X, **params):
        """Perform DBSCAN clustering from vector array or distance matrix.

        Parameters
        ----------
        X: array [n_samples, n_samples] or [n_samples, n_features]
            Array of distances between samples, or a feature array.
            The array is treated as a feature array unless the metric is
            given as 'precomputed'.
        """
        if (self.indexing == None):
            # automatically choose based on photon count.  Use indexing automatically for high photons counts
            if len(X[0])>2000:
                self.core_sample_indices_, self.labels_ = dbscan3(X,
                                                         **self.get_params())
            else:
                self.core_sample_indices_, self.labels_ = dbscan(X,
                                                         **self.get_params())
        if (self.indexing ==True):
            self.core_sample_indices_, self.labels_ = dbscan3(X,
                                                         **self.get_params())
        elif (self.indexing == False):
            self.core_sample_indices_, self.labels_ = dbscan(X,
                                                         **self.get_params())
        return self

    