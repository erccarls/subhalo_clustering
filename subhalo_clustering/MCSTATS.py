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


def DBSCAN_Compute_Clusters(mcSims, eps = 0.1875,min_samples = 3,nCorePoints = 3, S_cut=2.0,numAnalyze=0, fileout = '',BGTemplate = 'BGRateMap.pickle', HESS = False, angularSize = 10.,SNR= 0.75,numProcs = 10, indexing = True):
    '''
    Main DBSCAN cluster method.  Input a list of simulation outputs and output a list of clustering properties for each simulation.
    Inputs:
        mcSims: this is the output from the runMC methods (i.e. the list of simulation outputs) (python object not filename)
        eps: DBSCAN epsilon parameter
        min_samples: min num points in epislon neighborhood for the DBSCAN algorithm.
        nCorePoints: After DBSCAN is run, there must be at least this many points for a cluster to not be thrown out.
        S_cut: Clusters used in statistics must be at least this significance level
        numAnalyze: number of simulations to analyze out of list.  default is 0 which analyzes all of them
        fileout: if not empty string, determines the file to store all the clustering info.
        BGTemplate: background template filename.
        -SNR is the percentage of photons that are background
        
    return: (cluster_Scale, cluster_S, cluster_Count, cluster_Members, cluster_stdevs) as described in draft 
        -cluster_stdevs is the significance weighted RMS of the clustering scale 

    '''
    
    # Initialize the thread pool
    if (numProcs<=0):numProcs += mp.cpu_count()
    p = pool.Pool(numProcs)
    
    # Check number to analyze
    if ((numAnalyze == 0) or (numAnalyze > len(mcSims))):
        numAnalyze =len(mcSims)
    print 'Analyzing ' + str(numAnalyze) + ' simulations using ' , numProcs, " CPUs..."
    
    # Define a single input function callable by each thread (async map can only take one argument)
    def DBSCAN_THREAD(sim, eps, min_samples,nCorePoints,indexing= None):
        X = zip(sim[0],sim[1])
        return (DBSCAN.RunDBScan(X, eps, min_samples,nCorePoints = nCorePoints, indexing = indexing), len(X[0]))
    DBSCAN_PARTIAL = partial(DBSCAN_THREAD,  eps=eps, min_samples=min_samples,nCorePoints = nCorePoints)
    
    # Call async map. 
    dbscanResults = p.map(DBSCAN_PARTIAL, mcSims[:numAnalyze])

    # Serial Version.  Only use for debugging
    #dbscanResults = map(DBSCAN_PARTIAL, mcSims[:numAnalyze])
    
    # Single Call Version. Useful for Debugging
    #dbscanResults = DBSCAN_THREAD(mcSims[0],  eps=eps, min_samples=min_samples,nCorePoints = nCorePoints)

    
    BGTemplate = pickle.load(open(BGTemplate,'r'))

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
        S = [DBSCAN.Compute_Cluster_Significance(cluster, BGTemplate, numPhotons,angularSize = angularSize,SNR = SNR) for cluster in clusters]
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
        #cluster_S.append(np.average(clusterSigs,weights = clusterMembers))
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
    


def Compute_Radial_Profiles(MCSims,numBins = 20,rmax=5.0,label = ''):
    """
    Computes radial event density profiles plots (from center of simulation)
    
    Input:
        -MCSims: simulations list
        -bins: radial bins, linearly spaced
        -rmax: maximum radius to use. beware of edge effects.
    """
    
    
    ppa = 300./10. # pixels per degree
    
    bins = np.linspace(0, 5, numBins+1)
    bgtemplate = pickle.load(open('BGRateMap.pickle','r'))
    
    BG = []
    #print DBSCAN.Evaluate_BG_Contribution(0.0, 0.0, 10*ppa, BGTemplate=bgtemplate, numBGEvents=.75*500, flatLevel=0)
    hists = []
    for i in range(numBins):
        hists.append([])
        #bgevents = DBSCAN.Evaluate_BG_Contribution(0.0, 0.0, bins[i+1]*ppa, BGTemplate=bgtemplate, numBGEvents=.75*500, flatLevel=0)-DBSCAN.Evaluate_BG_Contribution(0.0, 0.0, bins[i]*ppa, BGTemplate=bgtemplate, numBGEvents=.75*500, flatLevel=0)
        #BG.append(bgevents)
        #print BG 
    for sim in MCSims:
        r_list = []
        
        xVec = sim[0]
        yVec = sim[1]
        
        
        
        for i in range(len(xVec)):
            # Find radius from center
            r_list.append(math.sqrt(xVec[i]**2.0+yVec[i]**2.0))
        
        values = np.histogram(r_list, bins=bins)[0]
        
        # Need to compute density
        for i in range(numBins):
            area = math.pi* (bins[i+1]**2.0-bins[i]**2.0)
            values[i] = float(values[i])/area
            
        for i in range(20):
            hists[i].append(values[i])
        
    Y, YERR = [],[]
    for i in range(20):
        Y.append(np.mean(hists[i]))
        YERR.append(np.std(hists[i]))
    
    
    #plt.step(bins[:-1], values)
    plt.errorbar(bins[:-1], Y, YERR,marker='o',label = label)
    plt.yscale('log')
    plt.xscale('log')
    
    return bins[:-1], Y, YERR
            
            
            
#===============================================================================
# More plotting tools that are older.
#===============================================================================
def Plot_Rate_Map(mapName,angularSize,fileOut):
    """Plot the map of annihilation rate."""
    
    map = pickle.load(open(mapName, "r" ))
    img = plt.imshow(map,origin='lower', extent=[-angularSize/2,angularSize/2,-angularSize/2,angularSize/2])
    plt.colorbar(orientation='vertical')
    
    plt.xlabel(r'$l[^\circ]$')
    plt.ylabel(r'$b[^\circ]$')
    plt.title(r'$\rho_{DM}^2$')

    plt.savefig(str(fileOut)+ '.png')
    plt.show()
    plt.clf()
    
    
def Plot_MC_Positions(MC,fileOut='',angularSize = 10.0):
    """
    Plot results of a Monte Carlo simulation.
        input:
        -MC: One Monte Carlo simulation.  If taking from a set of runs, pass mcSims[i] for the i'th simulation
        -fileOut: if not blank string then saves the file to this name
    """
    plt.clf()
    plt.scatter(MC[0], MC[1],color = 'b',s=4,marker='+')
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

    plt.show()    


def Plot_Fermi_Data(x,y,fileOut='',xLim=(-5,5),yLim=(-5,5), colormap=False, t = []):
    """
    Scatter plot of points
        input: 
            x: vector of x coords
            y: vector of y coords
        can also plot time dependence, but really no point.
    """
    
    plt.clf()
    fig = plt.figure(1)
    ax1 = fig.add_subplot(1,1,1, aspect='equal')
    
    if (colormap == True and t != []):
        scatter = ax1.scatter(x,y,c=t)
        plt.colorbar(scatter)
    else:
        ax1.scatter(x, y)#, s=5,lw =2,color = 'b')
    ax1.set_xlim(xLim[1],xLim[0])
    ax1.set_ylim(yLim[0],yLim[1])
    
    plt.xlabel(r'$l[^\circ]$')
    plt.ylabel(r'$b[^\circ]$')
    
    if fileOut != '':
        pp = PdfPages(fileOut + '.pdf')
        plt.savefig(pp, format='pdf')
        print "Figures saved to ", str(fileOut)+ '.pdf\n',
        pp.close()
    plt.show()    

        
def Plot_MC_Centroids(X,Y,fileOut,xLim=(2.5,-2.5),yLim=(-2.5,2.5)):
    """Plot centroids from a list of monte carlo simulations.  Needs updating"""
    import pylab,math
    plt.clf()
    
    fig = plt.figure(1)
    fig.add_subplot(1,1,1,aspect = 'equal')
    plt.scatter(X,Y, s=1,lw =0,color = 'b',rasterized=True)
    
    plt.xlabel(r'$l[^\circ]$')
    plt.ylabel(r'$b[^\circ]$')
    #plt.title(r'Centroids')
    
    r = []
    for i in range(len(X)):
        r.append(math.sqrt(X[i]**2.+Y[i]**2.))
    sigma = np.std(r)
    cir = pylab.Circle((np.average(X),np.average(Y)), radius=sigma, alpha =.3, fc='b')
    cir2 = pylab.Circle((np.average(X),np.average(Y)), radius=sigma*2., alpha =.3, fc='r')
    cir3 = pylab.Circle((np.average(X),np.average(Y)), radius=sigma*3., alpha =.3, fc='y')
    
    cent1 = (-0.12672151515151459, 0.32554675757575757) #106-121
    cent2 = (-0.768305075000001, 0.13957773025000009) # 121-136
    cent3 = (-1.2090912941176524, 0.28272128823529419) #136-151
    plt.scatter(cent1[0],cent1[1], marker = '+', s = 200,lw=2, color = 'r')
    plt.scatter(cent2[0],cent2[1], marker = '+', s = 200,lw=2, color = 'g')
    plt.scatter(cent3[0],cent3[1], marker = '+', s = 200,lw=2, color = 'purple')
    
    plt.legend(('MC','106-121 GeV','121-136 GeV','136-151 GeV'))
    pylab.gca().add_patch(cir3)
    pylab.gca().add_patch(cir2)
    pylab.gca().add_patch(cir)
    
    plt.xlim(xLim)
    plt.ylim(yLim)
    
    if fileOut!='':
        pp = PdfPages(fileOut + '.pdf')
        plt.savefig(pp, format='pdf',dpi=300)
        print "Figures saved to ", str(fileOut)+ '.pdf\n',
        pp.close()

    #plt.show()    
    
    


    


        
            
            
            
    
    
    
    
    
    