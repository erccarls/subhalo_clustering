'''
Statistics and plotting tools for Fermi MC tools

Created on Jul 24, 2012
@author: Eric Carlson
'''

import matplotlib.pyplot as plt  #@UnresolvedImport
import matplotlib.image as mpimg #@UnresolvedImport
import matplotlib.cm as cm #@UnresolvedImport
import matplotlib, scipy #@UnresolvedImport
import pickle, sys
import numpy as np
import scipy.cluster as cluster
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.mlab as mlab
import math, time
from scipy.cluster.hierarchy import *
from operator import itemgetter
import DBSCAN

import multiprocessing as mp
from multiprocessing import pool
from functools import partial


def DBSCAN_Compute_Clusters(mcSims, eps, min_samples ,nCorePoints = 3, numAnalyze=0, fileout = '',numProcs = 1):
    '''
    Main DBSCAN cluster method.  Input a list of simulation outputs and output a list of clustering properties for each simulation.
    Inputs:
        mcSims: this is the output from the runMC methods (i.e. the list of simulation outputs) (python object not filename)
        eps: DBSCAN epsilon parameter
        min_samples: min num points in epislon neighborhood for the DBSCAN algorithm.
    Optional Inputs:
        nCorePoints=3: After DBSCAN is run, there must be at least this many points for a cluster to not be thrown out.
        numAnalyze=0 : number of simulations to analyze out of list.  default is 0 which analyzes all of them
        fileout=''   : if not empty string, store all the clustering info in a pickle file.
    Returns:
        dbScanResults: a tuple (clusterReturn, labels) for each simulation
            clusterReturn: For each cluster, a list of points in that cluster
            labels: A list of points with cluster identification number.  -1 corresponds to noise.
                    NOTE: ALL BORDER POINTS ARE CONSIDERED NOISE.  See DBSCAN.py for info if you need
                    Border points.
    '''
    
    # Initialize the thread pool
    if (numProcs<=0):numProcs += mp.cpu_count()
    p = pool.Pool(numProcs)
    
    # Check number to analyze
    if ((numAnalyze == 0) or (numAnalyze > len(mcSims))):
        numAnalyze =len(mcSims)
    print 'Analyzing ' + str(numAnalyze) + ' simulations using ' , numProcs, " CPUs..."
    
    

    DBSCAN_PARTIAL = partial(DBSCAN_THREAD,  eps=eps, min_samples=min_samples,nCorePoints = nCorePoints)
    
    # Call mutithreaded map. 
    dbscanResults = p.map(DBSCAN_PARTIAL, mcSims[:numAnalyze])

    # Serial Version.  Only use for debugging
    #dbscanResults = map(DBSCAN_PARTIAL, mcSims[:numAnalyze])
    
    # Single Call Version. Useful for Debugging
    #dbscanResults = DBSCAN_THREAD(mcSims[0],  eps=eps, min_samples=min_samples,nCorePoints = nCorePoints)
    
    # Write to file if requested
    if (fileout != ''): pickle.dump(dbscanResults, open(fileout,'wb'))
    
    return dbscanResults


# Define a single input function callable by each thread (async map can only take one argument)
def DBSCAN_THREAD(sim, eps, min_samples,nCorePoints,indexing= None):
    X = zip(sim[0],sim[1])
    return DBSCAN.RunDBScan(X, eps, min_samples,nCorePoints = nCorePoints, indexing = indexing)



def Cluster_Sigs_BG(dbscanResults, BGTemplate = 'BGRateMap.pickle',angularSize = 10.,BG= 0.75,numProcs = 1):
    """
    Compute the cluster significances on results of DBSCAN_Compute_Clusters() using a background model.
    
    Inputs:
        dbscanResults: output from DBSCAN_Compute_Clusters.  Must load from file if using pickled results
        BGTemplate: background template filename.
        angularSize: Size of square in degrees
        BG: The expected percentage of photons that are background
    """
    # Initialize the thread pool to the correct number of threads
    if (numProcs<=0):numProcs += mp.cpu_count()
    p = mp.pool.Pool(numProcs)
    
    # Load background template
    BGTemplate = pickle.load(open(BGTemplate,'r'))

    # Asynchronosly map results
    BG_PARTIAL = partial(BG_THREAD, BGTemplate= BGTemplate, angularSize = angularSize, BG = BG)
    return p.map(BG_PARTIAL,dbscanResults)

def BG_THREAD(sim,BGTemplate, angularSize , BG ):
    clusters,labels = sim 
    return [DBSCAN.Compute_Cluster_Significance(cluster, BGTemplate, len(labels),angularSize = angularSize,BG = BG) for cluster in clusters]


#===============================================================================
# DEPRECATED 
#===============================================================================
def Profile_Clusters(dbscanResults, BGTemplate = 'BGRateMap.pickle',S_cut=2.0,angularSize = 10.,BG= 0.75, fileout=''):
    """
    Computes properties on the results of DBSCAN_Compute_Clusters()
        input:
            dbscanResults: output from DBSCAN_Compute_Clusters.  Must load from file if using pickled results
        Optional Inputs
            S_cut: Clusters used in statistics must be at least this significance level
            BGTemplate = String with path to the background template file
            BG: The expected percentage of photons that are background
            
        
        Returns: (cluster_Scale, cluster_S, cluster_Count, cluster_Members, cluster_stdevs) as described in draft 
            cluster_stdevs is the significance weighted RMS of the clustering scale 
    """
    
    
    cluster_S = []       # Mean Significance weighted by number of cluster members for ALL CLUSTERS    
    cluster_Count = []   # Mean Number of clusters found s>s_cut
    cluster_Scale = []   # Mean Cluster Scale weighted by significance for S>s_cut
    cluster_Members = [] # Mean number of cluster Members s> s_cut
    cluster_Out = []     # for each sim a tuple of (cluster_s, num_members) for each cluster in the simulation.
     
    for sim in dbscanResults: 
        clusters,labels = sim
        numPhotons = len(labels) 
        #===========================================================================
        # Compute Cluster Properties
        #===========================================================================     
        # Determine Cluster Significance from background template
        S = [DBSCAN.Compute_Cluster_Significance(cluster, BGTemplate, numPhotons,angularSize = angularSize,SNR = BG) for cluster in clusters]
        # Number of clusters
        clusterMembersAll = [len(cluster) for cluster in clusters]
        # list of pairs (s,num members) for each cluster
        cluster_Out.append(zip(S, clusterMembersAll))
        # S>S_Cut Cluster Indexes
        sigClustersIDX = np.where((S>=S_cut))[0]
        # S>S_Cut Clusters
        sigClusters = clusters[clusters]
        
        sigs = S[sigClustersIDX]
        # Compute Cluster Scales
        scale = [DBSCAN.Compute_Cluster_Scale(cluster)[0] for cluster in sigClusters]
        # S>S_Cut Cluster Member Counts
        members = clusterMembersAll[sigClustersIDX]
        
        
        #===========================================================================
        # Compute Weighted Means and append to master list 
        #===========================================================================
        # Append All cluster sigs.  Rest of quantities require S>2.0
        cluster_S.append(np.average(S,weights = clusterMembersAll))
        cluster_Count.append(len(sigs))
        if len(sigs)!=0:
            cluster_Scale.append(np.average(scale, weights = sigs))
            cluster_Members.append(np.average(members, weights = sigs))
        
    output = (cluster_Scale, cluster_S, cluster_Count, cluster_Members, cluster_Out)
    # Write results to file
    if fileout != '':
        pickle.dump(output, open(fileout, 'wb'))
    return output



#################################################################################################
# Plotting Tools
#################################################################################################

#TODO: 
def Plot_Cluster_Scales(models, labels,xlabel, fig, subplot, bins = 100, fileout = '', PlotFermiScale = False):
    """
    Generates the results summary plots.  See sample usage in runMC_v2
    """
    
    width = .08
      
    fig.add_subplot(4,2,abs(subplot))      
    
    hist = []
    for i in models:
        if (abs(subplot) == 5 or abs(subplot) == 6):
            if subplot>0:
                hist.append(np.histogram(i, bins=np.linspace(0, 9, 10)))
            elif subplot<0:
                hist.append(np.histogram(i, bins=np.linspace(0, 9, 10)))
        else:
            hist.append(np.histogram(i, bins=bins))
        
    
    c = ['b','g','r','c','m','y','k','b','g','r']
    
    for i in range(len(hist)):
        if  (abs(subplot) == 6 or abs(subplot) == 5):
            #plt.step(hist[i][1][:-1], np.array(hist[i][0],'float')/len(models[i]), label=labels[i])
            if subplot>0 and 'Pulsar' in labels[i]:
                plt.bar(hist[i][1][:-1]+width*(i-7), np.array(hist[i][0],'float')/len(models[i]) , width ,fill = True, label=labels[i], color =c[i])
            elif 'Pulsar' in labels[i]:
                plt.bar(hist[i][1][:-1]+width*(i-7), np.array(hist[i][0],'float')/len(models[i]) , width ,fill = True, color =c[i])
            else:
                plt.bar(hist[i][1][:-1]+width*(i-7), np.array(hist[i][0],'float')/len(models[i]) , width ,fill = False, label=labels[i], edgecolor =c[i])
            
            
        elif 'Pulsar' in labels[i] and subplot > 0:
            plt.step(hist[i][1][:-1], np.array(hist[i][0],'float')/len(models[i]), label=labels[i])
        #elif 'Pulsar' in labels[i] and subplot < 0:
        #    plt.step(hist[i][1][:-1], np.array(hist[i][0],'float')/len(models[i]), label=labels[i], ls = '--')
        elif abs(subplot) != 5 and abs(subplot) != 6:
            plt.step(hist[i][1][:-1], np.array(hist[i][0],'float')/len(models[i]),ls ='--', label=labels[i])
        else:
            continue
            
#        elif abs(subplot) != 5 and abs(subplot) != 6: 
#            plt.step(hist[i][1][:-1], np.array(hist[i][0],'float')/len(models[i]),ls ='--', label=labels[i])
#        else: 
#            return
    if subplot == 1 or subplot == 2:
        high = plt.ylim()[1]
        ############################################
        # For Sideband bkg model
#        plt.fill_betweenx([0,10], [1.198,1.198], [1.44,1.44], facecolor='m', alpha=0.05)
#        plt.axvline(1.337,ls='-.',color = 'm',label = 'Fermi 120-140 GeV')
        ############################################
        # For Galprop bkg model
        #plt.fill_betweenx([0,10], [1.287,1.287], [1.56,1.56], facecolor='m', alpha=0.05)
        #plt.axvline(1.4429,ls='-.',color = 'm',label = 'Fermi 120-140 GeV')
        if subplot == 1:
            plt.axvline(1.287,ls='-.',color = 'm',label = 'Fermi 120-140 GeV')
        
        plt.ylim(0,.25)
        if subplot == 5:
            plt.xlim(0,14)
        
    if  subplot == 3 and PlotFermiScale ==True:
        high = plt.ylim()[1]
        #plt.fill_betweenx([0,10], [.2192,.2192], [.439,.439], facecolor='m', alpha=0.05)
        ############################################
        # For Sideband bkg model
        #plt.axvline(.339,ls='-.',color = 'm')
        ############################################
        # For Galprop bkg model
        #plt.axvline(.34,ls='-.',color = 'm')
        plt.axvline(.2192,ls='-.',color = 'm')
        
        plt.ylim(0,.3)
        plt.xlim(0,.5)
        
    if subplot == 4:
        plt.ylim(0,.3)
        plt.xlim(0,.06)
    
        
    if (subplot == 5) and  PlotFermiScale ==True:
        #plt.axvline(2,ls='-.',color = 'm')
        plt.xlim(-1,8)
        plt.ylim(0,1)
        
        
    if subplot==6:
        plt.xlim(-1,8)
        plt.ylim(.01,1)
        #plt.yscale('log')
        
    
    if subplot == 7 and PlotFermiScale ==True:
        plt.axvline(4,ls='-.',color = 'm')
        plt.ylim(0,.3)
    if subplot == 8:
        plt.ylim(0,.3)
        
    plt.xlabel(xlabel)
    if subplot in [1,3,5,7]:
        plt.ylabel(r'$f$')
    



            
            
            
#===============================================================================
# More plotting tools that are older.
#===============================================================================
def Plot_Rate_Map(mapName,angularSize,fileOut):
    """
    Plot the map of annihilation rate.
        inputs:
            mapName: string with filename to pickle file of rate map.
            
    
    """
    plt.clf()
    map = pickle.load(open(mapName, "r" ))
    img = plt.imshow(map,origin='lower', extent=[-angularSize/2,angularSize/2,-angularSize/2,angularSize/2])
    plt.colorbar(orientation='vertical')
    
    plt.xlabel(r'$l[^\circ]$')
    plt.ylabel(r'$b[^\circ]$')
    plt.title(r'$\rho_{DM}^2$')

    plt.savefig(str(fileOut)+ '.png')
    plt.show()
    

#===============================================================================
# Quick way to plot the monte-carlo or Fermi data, just pass a tuple of coordinate
# pairs for each positions. 
#===============================================================================
def Plot_MC_Positions(MC,fileOut='',angularSize = None):
    """
    Plot results of a Monte Carlo simulation.
        input:
        -MC: One Monte Carlo simulation.  If taking from a set of runs, pass mcSims[i] for the i'th simulation
        -fileOut: if not blank string then saves the file to this name and does not display plot
        -angularSize: If None, auto-scales to the min/max of simulation.  Otherwise pass the width in simulation coordinates
    """
    plt.clf()
    plt.scatter(MC[0], MC[1],color = 'b',s=4,marker='+')
    if angularSize != None:
        plt.xlim(angularSize/2.0,-angularSize/2.0)
        plt.ylim(-angularSize/2.0,angularSize/2.0)
            
        plt.xlabel(r'$l[^\circ]$')
        plt.ylabel(r'$b[^\circ]$')
        plt.title(r'Count Map')
        if fileOut!='':
            pp = PdfPages(fileOut + '.pdf')
            plt.savefig(pp, format='pdf')
            print "Figures saved to ", str(fileOut)+ '.pdf\n',
            pp.close()
            plt.savefig(fileOut + '.png', format='png')
            return
        plt.show()    

