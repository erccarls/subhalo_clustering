import pickle
import MC,MCSTATS,DBSCAN
import numpy as np
import sys, math
from matplotlib.backends.backend_pdf import PdfPages

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

#angularSize = 2 # box surrounding galactic center
#size = 200  # size of output ratemap
#profile = ('EIN',20.0,0.17)
#fileOut = 'ein_test.pickle'
#map = MC.Gen_Annihilation_Map(angularSize, size, profile,fileOut)




def run_MC_NFW():
    """Generate Rate Map and Run Monte-Carlo NFW""" 
    #NFW Rate Map
    profile = ('NFW',23.5,1.0)
    if run48:
        fileOut = 'NFWRateMap.pickle'
        #map = MC.Gen_Annihilation_Map(angularSize, size, profile,fileOut)
        #MCSTATS.Plot_Rate_Map('NFWRateMap.pickle', angularSize, fileOut)
        # MC
        mcSims = MC.RUN(numSims, fileOut, numPhotons=numPhotons,angularSize=angularSize, outputSize=outputSize, mcList='MCOut_galprop_NFW_48.pickle',numProcs = numProcs)
    
    if run2000:
        fileOut = 'NFWRateMapHESS.pickle'
        #map = MC.Gen_Annihilation_Map(angularSize2, size, profile,fileOut)
        mcSims = MC.RUN(numSims2, fileOut, numPhotons=numPhotons2,angularSize=angularSize2, outputSize=outputSize, mcList='MCOut_galprop_NFW_500.pickle',HESS=True,numProcs = numProcs)


def run_MC_EIN():
    """Generate Rate Map and Run Monte-Carlo Einasto"""
    # Einasto rate map
    profile = ('EIN',20,0.17)
    if run48:
        fileOut = 'EINRateMap.pickle'
        #map = MC.Gen_Annihilation_Map(angularSize, size, profile,fileOut)
        #MCSTATS.Plot_Rate_Map('EINRateMap.pickle', angularSize, fileOut)
        # MC
        mcSims = MC.RUN(numSims, fileOut, numPhotons=numPhotons,angularSize=angularSize, outputSize=outputSize, mcList='MCOut_galprop_EIN_48.pickle',numProcs = numProcs)
    if run2000:
        fileOut = 'EINRateMapHESS.pickle'
    #    #map = MC.Gen_Annihilation_Map(angularSize2, size, profile,fileOut)
        mcSims = MC.RUN(numSims2, fileOut, numPhotons=numPhotons2,angularSize=angularSize2, outputSize=outputSize, mcList='MCOut_galprop_EIN_500.pickle',HESS=True,numProcs = numProcs)
    
def run_MC_FLAT():
    """Generate Rate Map and Run Monte-Carlo Flat"""
    # Einasto rate map
    profile = ('FLAT',0,0)
    if run48:
        fileOut = 'FLATRateMap.pickle'
        #map = MC.Gen_Annihilation_Map(angularSize, size, profile,fileOut)
        #MCSTATS.Plot_Rate_Map('FLATRateMap.pickle', angularSize, fileOut)
        # MC
        mcSims = MC.RUN(numSims, fileOut, numPhotons=numPhotons,angularSize=angularSize, outputSize=outputSize, mcList='MCOut_galprop_FLAT_48.pickle',numProcs = numProcs)
    if run2000:
        fileOut = 'FLATRateMapHESS.pickle'
        #map = MC.Gen_Annihilation_Map(angularSize2, size, profile,fileOut)
        mcSims = MC.RUN(numSims2, fileOut, numPhotons=numPhotons2,angularSize=angularSize2, outputSize=outputSize, mcList='MCOut_galprop_FLAT_500.pickle',HESS=True,numProcs = numProcs)

def run_MC_NFWDECAY():
    """Generate Rate Map and Run Monte-Carlo NFW but not squared"""
    # NFW_DECAY rate map
    profile = ('NFWDECAY',23.5,1.0)
    if run48:
        fileOut = 'NFWDECAYRateMap.pickle'
        #map = MC.Gen_Annihilation_Map(angularSize, size, profile,fileOut)
        #MCSTATS.Plot_Rate_Map('NFWDECAYRateMap.pickle', angularSize, fileOut)
        # MC
        mcSims = MC.RUN(numSims, fileOut, numPhotons=numPhotons,angularSize=angularSize, outputSize=outputSize, mcList='MCOut_galprop_NFWDECAY_48.pickle',numProcs = numProcs)
    if run2000:
        fileOut = 'NFWDECAYRateMapHESS.pickle'
        #map = MC.Gen_Annihilation_Map(angularSize2, size, profile,fileOut)
        mcSims = MC.RUN(numSims2, fileOut, numPhotons=numPhotons2,angularSize=angularSize2, outputSize=outputSize, mcList='MCOut_galprop_NFWDECAY_500.pickle',HESS=True,numProcs = numProcs)
    
def run_MC_PULSAR(numPulsars=6):
    """
    Generate Rate Map and Run Monte-Carlo NFW but not squared
        -input: 
            numPulsars: number 
    """
    # Pulsars
    profile = ('PULSAR',)
    if run48:
        fileOut = 'PULSARRateMap.pickle'
        #map = MC.Gen_Annihilation_Map(angularSize, size, profile,fileOut)
        #MCSTATS.Plot_Rate_Map('PULSARRateMap.pickle', angularSize, fileOut)
        # MC
        mcSims = MC.RUN_PULSAR(numSims, fileOut, numPhotons=numPhotons,numPulsars=numPulsars,angularSize=angularSize2, outputSize=outputSize, mcList='MCOut_galprop_PULSAR_48_' + str(numPulsars) + '.pickle',numProcs = numProcs)
    if run2000:
        fileOut = 'PULSARRateMapHESS.pickle'
        #map = MC.Gen_Annihilation_Map(angularSize2, size, profile,fileOut)
        mcSims = MC.RUN_PULSAR(numSims2, fileOut, numPhotons=numPhotons2,numPulsars=numPulsars,angularSize=angularSize2, outputSize=outputSize, mcList='MCOut_galprop_PULSAR_500_' + str(numPulsars) + '.pickle',HESS=True,numProcs = numProcs)
    
    
def run_DBSCAN_ANALYSIS(mcFile,fileout, numAnalyze = 0, eps=.201719, S_cut = 2.0,nCore = 3, nCluster = 3,angularSize = angularSize, HESS = False,SNR = .75):
    """Run cluster analysis on mcFile output from MC.RUN."""
    print 'Loading MC File...'
    mcSims = pickle.load(open(mcFile, "r" ))
    print 'Running Cluster Analysis...'
    if HESS ==True:
        clusterResults = MCSTATS.DBSCAN_Compute_Clusters(mcSims, S_cut=S_cut, eps=eps , numAnalyze = numAnalyze, fileout = fileout,n_cluster=nCluster, nCore = nCore,angularSize = angularSize2,HESS = HESS,SNR= SNR,numProcs = numProcs)
    else:
        clusterResults = MCSTATS.DBSCAN_Compute_Clusters(mcSims, S_cut=S_cut, eps=eps , numAnalyze = numAnalyze, fileout = fileout,n_cluster=nCluster, nCore = nCore,angularSize = angularSize,HESS = HESS,SNR= SNR,numProcs = numProcs)
    print "File Output to " + fileout + '...'
    return

#run_MC_PULSAR(numPulsars=6)

# Testing Code
#mcSims = MC.RUN(50, 'NFWRateMap.pickle', numPhotons=500,angularSize=angularSize, outputSize=outputSize, mcList='test.pickle' + str(1) + '.pickle',flatLevel = 0.0)
#MCSTATS.Compute_Radial_Profiles(mcSims, bins=20)

#print MCSTATS.DBSCAN_Compute_Clusters(mcSims, flatLevel=0.00, numAnalyze = 0)


# --------------------- Test Plotting
profile = ('PULSAR',)
fileOut = 'PULSARRateMap.pickle'
#mcSims = MC.RUN_PULSAR(1000, fileOut, numPhotons=numPhotons,numPulsars=4,angularSize=angularSize, outputSize=outputSize, mcList='test.pickle',numProcs=6)
##print MCSTATS.DBSCAN_Compute_Clusters(mcSims, S_cut=2.0, eps=.15 , numAnalyze = 0, fileout = '',n_cluster=3, nCore = 3,angularSize = 10.0)
##MCSTATS.Plot_MC_Positions(mcSims[0],angularSize = 10)
#
#results = MCSTATS.DBSCAN_Compute_Clusters(mcSims, eps=.15, n_cluster=3, nCore=3, S_cut=2.0, numAnalyze=0,  HESS=False, angularSize=angularSize,numProcs = 6)
#count = 0
#
#clus = []
#for sim in results[4]:
#    for cluster in sim:
#        clus.append(cluster)
#        if cluster[0] > 1.287: 
#            count+=1 
#            break
#print count
        


fileOut = 'NFWRateMapHESS.pickle'
profile = ('NFW',23.5,1.0)
#map = MC.Gen_Annihilation_Map(angularSize2, size, profile,fileOut)
#mcSims = MC.RUN(20, fileOut, numPhotons=numPhotons2,angularSize=angularSize2, outputSize=outputSize, mcList='MCOut_galprop_NFW_500.pickle',HESS=True)
#MCSTATS.Plot_MC_Positions(mcSims[0], 'test', angularSize2)

    

profile = ('PULSAR',)
fileOut = 'PULSARRateMapHESS.pickle'
#map = MC.Gen_Annihilation_Map(angularSize2, size, profile,fileOut)
#mcSims = MC.RUN_PULSAR(numSims, fileOut, numPhotons=100,numPulsars=numPulsars,angularSize=angularSize, outputSize=outputSize, mcList='MCOut_galprop_PULSAR_100_' + str(numPulsars) + '.pickle')
#fileOut = 'PULSARRateMapHESS.pickle'
#mcSims = MC.RUN_PULSAR(1, fileOut, numPhotons=2000,numPulsars=1,angularSize=angularSize2, outputSize=outputSize, mcList='test.pickle',HESS=True)
#MCSTATS.Plot_MC_Positions(mcSims[0], 'test', angularSize = angularSize2)
#MCSTATS.DBSCAN_Compute_Clusters(mcSims, eps=.16, n_cluster=3, nCore=5, S_cut=2.0, numAnalyze=0,  HESS=True, angularSize=angularSize2)


import matplotlib.pyplot as plt 
fig2 = plt.figure(2, figsize = [12,8])
fig = plt.figure(1, (12,16))


def Gen_Plots2(fileout = 'Cluster_All_48'):
    """48 photon plots"""
    inputList = (('NFW','Clusters_MCOut_galprop_NFW_48.pickle'),
                 ('Flat','Clusters_MCOut_galprop_FLAT_48.pickle'), 
                 ('Ein','Clusters_MCOut_galprop_EIN_48.pickle'),
                 ('NFW Decay','Clusters_MCOut_galprop_NFWDECAY_48.pickle'),
                 ('Pulsar 1 ','Clusters_MCOut_galprop_PULSAR_48_1.pickle'),
                 ('Pulsar 2','Clusters_MCOut_galprop_PULSAR_48_2.pickle'),
                 ('Pulsar 3','Clusters_MCOut_galprop_PULSAR_48_3.pickle'),
                 ('Pulsar 4','Clusters_MCOut_galprop_PULSAR_48_4.pickle'),
                 ('Pulsar 5','Clusters_MCOut_galprop_PULSAR_48_5.pickle'),
                 ('Pulsar 6','Clusters_MCOut_galprop_PULSAR_48_6.pickle')
                 )
    labels = []
    scales = []
    sigs   = []
    count  = []
    members= []
    countLess = []
    
    for i in inputList:
        (cluster_Scale, cluster_S, cluster_Count, cluster_Members,cluster_Out) = pickle.load(open(i[1],'r'))
        scales.append(cluster_Scale)
        sigs.append(cluster_S)
        #count.append(cluster_Count)
        members.append(cluster_Members)
        labels.append(i[0])
        
        
        num_less= 0 # number of simulations with at least one cluster s>2
        num_less2 = 0
        num_less3 = 0
        full_count = []
        num_less_s_list = []
        for sim in cluster_Out:
            
            count2 = 0
            full_count.append(len(sim))
            num_less_s = 0    
            breakFlag = False
            for cluster in sim:
                if cluster[0] > 1.287:
                    num_less_s += 1
                if cluster[0] > 1.287 and breakFlag == False:
                    num_less += 1
                    breakFlag = True
                if cluster[0] > 1.287:
                    count2 += 1
                    if count2 == 2:
                        num_less2 += 1
                    if count2 == 3:
                        num_less3 += 1
                
            
            num_less_s_list.append(num_less_s)
        
        
        count.append(full_count)
        countLess.append(num_less_s_list)
        
        print i
        print 'Percentage of Sims with at least 1 cluster S > 1.287:' , float(num_less)/float(len(cluster_Out))
        print 'Percentage of Sims with at least 2 cluster S > 1.287:' , float(num_less2)/float(len(cluster_Out))
        print 'Percentage of Sims with at least 3 cluster S > 1.287:' , float(num_less3)/float(len(cluster_Out))
        print 
        
    plt.clf()
    
    #print len(sigs[5])
    
    from matplotlib.font_manager import FontProperties
    fontP = FontProperties()
    fontP.set_size('small')
    
    MCSTATS.Plot_Cluster_Scales(sigs, labels,r'S',fig,1,   bins=30, fileout='Cluster_Sigs',PlotFermiScale = True)
    
    
    MCSTATS.Plot_Cluster_Scales(scales, labels,r'$r_{cluster}[^\circ] | s > 1.287$',fig,3, bins=30, fileout='Cluster_Scaling',PlotFermiScale = True)
    
    MCSTATS.Plot_Cluster_Scales(countLess, labels,r'$N_{clusters} | s > 1.287$',fig,5,  bins=50, fileout='Cluster_Count',PlotFermiScale = True)
    #MCSTATS.Plot_Cluster_Scales(countLess, labels,r'$N_{clusters}$',fig,-5,  bins=50, fileout='Cluster_Count',PlotFermiScale = True)
    
    MCSTATS.Plot_Cluster_Scales(members, labels,r'$N_{members} | s > 1.287$',fig,7,bins=30, fileout='Cluster_Members',PlotFermiScale = True)
    
    if fileout != '':
        pp = PdfPages(fileout + '.pdf')
        plt.savefig(pp, format='pdf')
        print "Figures saved to ", str(fileout)+ '.pdf\n',
        pp.close()
    
    #plt.show()





def Gen_Plots3(fileout = 'Cluster_All_500'):
    """500 photon plots"""
    inputList = (('NFW','Clusters_MCOut_galprop_NFW_500.pickle'),
                 ('Flat','Clusters_MCOut_galprop_FLAT_500.pickle'), 
                 ('Ein','Clusters_MCOut_galprop_EIN_500.pickle'),
                 ('NFW Decay','Clusters_MCOut_galprop_NFWDECAY_500.pickle'),
                 ('Pulsar 1 ','Clusters_MCOut_galprop_PULSAR_500_1.pickle'),
                 ('Pulsar 2','Clusters_MCOut_galprop_PULSAR_500_2.pickle'),
                 ('Pulsar 3','Clusters_MCOut_galprop_PULSAR_500_3.pickle'),
                 ('Pulsar 4','Clusters_MCOut_galprop_PULSAR_500_4.pickle'),
                 ('Pulsar 5','Clusters_MCOut_galprop_PULSAR_500_5.pickle'),
                 ('Pulsar 6','Clusters_MCOut_galprop_PULSAR_500_6.pickle')
                 )
    labels = []
    scales = []
    sigs   = []
    count  = []
    members= []
    countLess = []
    
    for i in inputList:
        (cluster_Scale, cluster_S, cluster_Count, cluster_Members,cluster_Out) = pickle.load(open(i[1],'r'))
        scales.append(cluster_Scale)
        sigs.append(cluster_S)
        #count.append(cluster_Count)
        members.append(cluster_Members)
        labels.append(i[0])
        
        
        num_less= 0 # number of simulations with at least one cluster s>2
        num_less2 = 0
        num_less3 = 0
        full_count = []
        num_less_s_list = []
        for sim in cluster_Out:
            
            count2 = 0
            count3 = 0
            full_count.append(len(sim))
            num_less_s = 0    
            breakFlag = False
            for cluster in sim:
                if cluster[0] > 2.0:
                    num_less_s += 1
                if cluster[0] > 2.0 and breakFlag == False:
                    num_less += 1
                    breakFlag = True
                if cluster[0] > 2.0:
                    count2 += 1
                    if count2 == 2:
                        num_less2 += 1
                    if count2 == 3:
                        num_less3 += 1
                
            
            num_less_s_list.append(num_less_s)
        
        
        count.append(full_count)
        countLess.append(num_less_s_list)
        
        print i
        print 'Percentage of Sims with at least 1 cluster S > 2:' , float(num_less)/float(len(cluster_Out))
        print 'Percentage of Sims with at least 2 cluster S > 2:' , float(num_less2)/float(len(cluster_Out))
        print 'Percentage of Sims with at least 3 cluster S > 2:' , float(num_less3)/float(len(cluster_Out))
        print 
    
        
#        cluster20 = 0
#        cluster15 = 0 
#        for j in cluster_S:
#            if j<2.0:
#                cluster20 += 1.
#            if j<1.5:
#                cluster15 += 1.
        #print i, 'sig< 2.0', cluster20/len(cluster_S)  
        #print i, 'sig< 1.5', cluster15/len(cluster_S)
#        f =0.0
#        for j in cluster_Scale:
#            if j>.3397:
#                f+=1./float(len(cluster_Scale))
#        print i, ':',f
        
    #import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    fontP = FontProperties()
    fontP.set_size('small')
     
    #fig = plt.figure(2, (10,12))
    #plt.clf()
    MCSTATS.Plot_Cluster_Scales(sigs, labels,r'S',fig,2,   bins=30, fileout='Cluster_Sigs')
    #plt.legend(loc=1, ncol=1, fancybox=True, shadow=False,prop=fontP,borderaxespad=0.,labelspacing = .2)
    #plt.legend(loc=1, ncol=1, fancybox=True, shadow=False,prop=fontP,borderaxespad=0.,labelspacing = .2)
    MCSTATS.Plot_Cluster_Scales(scales, labels,r'$r_{cluster}[^\circ] | s > 2$',fig,4, bins=30, fileout='Cluster_Scaling')
    
    MCSTATS.Plot_Cluster_Scales(countLess, labels,r'$N_{clusters} | s > 2$',fig,6,  bins=30, fileout='Cluster_Count')
    #MCSTATS.Plot_Cluster_Scales(countLess, labels,r'$N_{clusters}$',fig,-6,  bins=30, fileout='Cluster_Count')
    
    MCSTATS.Plot_Cluster_Scales(members, labels,r'$N_{members} | s > 2$',fig,8,bins=30, fileout='Cluster_Members')
    plt.legend(loc=1, ncol=1, fancybox=True, shadow=False,prop=fontP,borderaxespad=0.,labelspacing = .2)
    
    if fileout != '':
        pp = PdfPages(fileout + '.pdf')
        plt.savefig(pp, format='pdf')
        print "Figures saved to ", str(fileout)+ '.pdf\n',
        pp.close()
    
    plt.show()
    #MCSTATS.Plot_Cluster_Contour(scales, sigs,labels, 'scale', 's')

def DBSCAN_PULSAR(numPulsars = 2, numPhotons =48,eps = 0.1875,n_cluster = 3,nCore = 3, S_cut=1.0, numSims = 100):
    
    BGTemplate = pickle.load(open('BGRateMap.pickle','r'))
    # Puslars
#    profile = ('PULSAR',)
#    fileOut = 'PULSARRateMap.pickle'
#    mcSims = MC.RUN_PULSAR(numSims, fileOut, numPhotons=numPhotons,numPulsars=numPulsars,angularSize=angularSize, outputSize=outputSize, mcList='test' + str(numPulsars) + '.pickle',flatLevel = 0.0,HESS = False)
    
    # NFW
    profile = ('NFW',23.5,1.0)
    fileOut = 'NFWRateMap.pickle'
    mcSims = MC.RUN(numSims, fileOut, numPhotons=numPhotons,angularSize=angularSize, outputSize=outputSize, mcList='test' + str(numPulsars) + '.pickle',HESS = False)
    
    cluster_Count = []   # Mean Number of clusters found
    cluster_Scale = []   # Mean Cluster Scale weighted by significance
    cluster_S = []       # Mean Significance weighted by number of cluster members
    cluster_Members = [] # Mean number of cluster Members
    
    
    
    for sim in mcSims:
        print sim[0]
        xVec = sim[0]
        yVec = sim[1]
        X = []
        for i in range(len(xVec)):
            X.append((xVec[i],yVec[i]))
        
        (clusters, clusterCount, noiseCount) = DBSCAN.RunDBScan(X, eps, n_cluster,nCore = nCore, plot=True)
        
    
        #===========================================================================
        # Compute Cluster Properties
        #===========================================================================
        clusterSigs = []
        clusterDist = []
        clusterSigmas = []
        clusterMembers = []
        for cluster in clusters:
            S =  DBSCAN.Compute_Cluster_Significance(cluster, BGTemplate, len(xVec),flatLevel = 0)
            if S>S_cut:
                d, sigma = DBSCAN.Compute_Cluster_Scale(cluster)  # compute the cluster Scale
                clusterSigs.append(S)
                clusterDist.append(d)
                clusterSigmas.append(sigma)
                clusterMembers.append(len(cluster))
                
                #print S, d, len(cluster)
        # If no clusters found
        if len(clusterSigs)==0:
            continue
        #===========================================================================
        # Compute Weighted Mean and Weighted Stdev 
        #===========================================================================
        clusterScale = np.average(clusterDist, weights = clusterSigs)
        weights = np.array(clusterSigs)/np.sum(clusterSigs)
        stdev = 0
        for i in range(len(clusterSigmas)):
            stdev += weights[i]**2 * clusterSigmas[i]**2 
        stdev = math.sqrt(stdev)
        #print clusterScale, stdev
        # Append to Master List
        cluster_Scale.append(clusterScale)
        cluster_S.append(np.average(clusterSigs,weights = clusterMembers))
        cluster_Count.append(len(clusters))
        cluster_Members.append(np.mean(clusterMembers))
        
    return cluster_Scale, cluster_S, cluster_Count, cluster_Members

#print DBSCAN_PULSAR(numSims = 1, numPhotons = 48,n_cluster = 3 , nCore =3,eps = .346)     
    
    
def DBSCAN_FERMI(energyRange = (120000,140000),eps=0.16,n_cluster = 3,nCore = 3, S_cut = 0.0):
    """Import Fermi Data and run DBSCAN"""
    import ParseFermi
    
    # Load fermi events 
    events = ParseFermi.Import_File('photons.txt', energyRange = energyRange,lonRange=(-5,5),latRange = (-5,5)) #@UndefinedVariable
    # Get into correct form
    X = []
    for i in events:
        X.append((i[1],i[2]))
        
    # Run DBSCAN
    (clusters, clusterCount, noiseCount) = DBSCAN.RunDBScan(X, eps, n_cluster,nCore = nCore,plot=False)
    # Load BG Template
    BGTemplate = pickle.load(open('BGRateMap.pickle','r'))
    
    #===========================================================================
    # Compute Cluster Properties
    #===========================================================================
    clusterSigsAll = []
    clusterMembersAll = []
    
    clusterSigs = []
    clusterDist = []
    clusterSigmas = []
    clusterMembers = []

    N_clust = 0
    for cluster in clusters:
        
        S =  DBSCAN.Compute_Cluster_Significance(cluster, BGTemplate, 48)
        
        # Keep a list of all the clusters found adn their weights
        clusterSigsAll.append(S)
        clusterMembersAll.append(len(cluster))
        
        print 'sig:'   , S
        print 'scale:' , DBSCAN.Compute_Cluster_Scale(cluster)
        print 'members', len(cluster)
        
        if S>S_cut:
            d, sigma = DBSCAN.Compute_Cluster_Scale(cluster)  # compute the cluster Scale
            clusterSigs.append(S)
            clusterDist.append(d)
            clusterSigmas.append(sigma)
            clusterMembers.append(len(cluster))
            #print S, d, len(cluster)
            
        if S>2.0:
            N_clust +=1
            
    # If no clusters found default to zero significance
    if len(clusters) == 0:
        clusterSigsAll.append(0.0)
        clusterMembersAll.append(1.0)
    
        
    return np.average(clusterDist), np.average(clusterSigsAll,weights = clusterMembersAll),len(clusters),np.mean(clusterMembers)


def epsPlot(subplot, filename, text):
    scale,S,count, members = pickle.load(open(filename,'rb'))
    #plt.subplot(2,3,i,adjustable='box', aspect=5)
    ax = plt.subplot(2,3,subplot)
    
    levels = []
    if subplot <=3:
        levels = np.linspace(0, 3, 13)
    else:
        levels = np.linspace(0, 8, 33)
    print levels
    
    CS = plt.contourf(S,levels = levels,extent = (EPSRANGE[0],EPSRANGE[-1],NMINRANGE[0], NMINRANGE[-1]),aspect='auto')
    plt.text(.08,3.5 ,text)
    
    plt.setp(ax.get_xticklabels(), rotation='vertical', fontsize=12)
    
    levels2 = list(np.linspace(1, 10, 10)) + list(np.linspace(12, 42, 11))
    CS2 = plt.contour(count,levels2,extent = (EPSRANGE[0],EPSRANGE[-1],NMINRANGE[0], NMINRANGE[-1]),colors = 'k')
    plt.clabel(CS2, inline=1, fontsize=12, colors = 'k')
    
    
    plt.ylabel(r'$N_{min}$')
    plt.xlabel(r'$\epsilon$')
    #if (subplot == 3 or subplot == 6):
    if ( subplot == 6):
        cax = fig2.add_axes([0.91, 0.1, 0.03, 0.365])
        cbar = plt.colorbar(CS , cax=cax)
        cbar.ax.set_ylabel('S')
    if ( subplot == 3):
        cax = fig2.add_axes([0.91, 0.535, 0.03, 0.365])
        cbar = plt.colorbar(CS , cax=cax)
        cbar.ax.set_ylabel('S')
    
    return
    



eps=0.1875         # 95% PSF containment for rear-converting events
eps=.15             # optimization

#print DBSCAN_FERMI( energyRange = (120000,140000),eps = eps, n_cluster = 3, nCore = 3)

#clusterData = DBSCAN_PULSAR(numPulsars = 3, numPhotons =48,eps = eps,n_cluster = 3, nCore = 2,numSims = 5, S_cut = 1.0)+ ('3 Pulsars',)
#MCSTATS.DBSCAN_STATS((clusterData,), fileOut='')
#run_DBSCAN_ANALYSIS('MCOut_galprop_PULSAR_48_2.pickle', 'clusters',numAnalyze =10)

if len(sys.argv) > 1: 
    print 'runMC.py beginning in mode ' + sys.argv[1]
    if (sys.argv[1] == '0'):
        run_MC_NFW()
    if (sys.argv[1] == '1'):
        run_MC_EIN()
    if (sys.argv[1] == '2'):
        run_MC_FLAT()
    if (sys.argv[1] == '3'):
        run_MC_NFWDECAY()
    if (sys.argv[1] == '4'):
        run_MC_PULSAR(int(sys.argv[2]))
        
    if (sys.argv[1] == '10'):
        run_DBSCAN_ANALYSIS(sys.argv[2],'Clusters_'+sys.argv[2], numAnalyze =0, S_cut = 1.287, nCore = 3, nCluster = 3,eps=.15)
    if (sys.argv[1] == '11'):
        # NOTE: nCore are required num core points for cluster.  NCluster is N_min required in epsilon neighb.
        SNR = .95
        run_DBSCAN_ANALYSIS(sys.argv[2],'Clusters_'+sys.argv[2], numAnalyze = 0, S_cut = 2.0, nCore = 3,nCluster = 8,eps=.25,angularSize = angularSize2, HESS = True,SNR = SNR)
        #Gen_Plots('Clusters_'+sys.argv[2])
    if (sys.argv[1] == '12'): 
        Gen_Plots2(fileout = 'Cluster_galprop_All_48')
        #Gen_Plots100(fileout = 'Cluster_galprop_All_100')
        Gen_Plots3(fileout = 'Cluster_galprop_All_500')



    #===========================================================================
    # SCAN SNR vs N_MIN
    #===========================================================================
    if (sys.argv[1] == '16'):
        profile = ('PULSAR',)
        fileOut = 'PULSARRateMap.pickle'
        
        numPulsars = 2
        numPhotons2 = 250
        angularSize = 10
        numSims = 20
        
        MCList = []
        for sig in np.arange(.0,.4,.025):
            MCList.append(MC.RUN_PULSAR(numSims, fileOut, numPhotons=numPhotons2,numPulsars=numPulsars,angularSize=angularSize, outputSize=outputSize, mcList='test_' + str(numPulsars) + '.pickle',HESS=True,Sig = sig))
        pickle.dump(MCList, open('MCList.pickle', 'wb'))
        
        
    if (sys.argv[1] == '17'):  
        """SCAN SNR vs N_MIN"""  
        MCList = pickle.load(open('MCList.pickle', 'rb'))
        numAnalyze = 20
        HESS = False
        
        
        NMINRANGE = range(3,20,1) 
        signal = np.arange(.0,.4,.025)
        
        eps = .1
        size= len(signal),len(NMINRANGE)
        clusterResults = [np.zeros(size),np.zeros(size),np.zeros(size),np.zeros(size)]
        
        for i in range(len(NMINRANGE)):
            for j in range(len(signal)):
                #print 'eps, nmin: ', signal[j], NMINRANGE[i]
                cluster_Scale, cluster_S, cluster_Count, cluster_Members, stdev = MCSTATS.DBSCAN_Compute_Clusters(MCList[j], S_cut=2.0, eps=eps , numAnalyze = numAnalyze,n_cluster=NMINRANGE[i], nCore = NMINRANGE[i],HESS=HESS)
                
                clusterResults[0][i][j] = np.mean(cluster_Scale)
                clusterResults[1][i][j] = np.mean(cluster_S)
                clusterResults[2][i][j] = np.mean(cluster_Count)
                clusterResults[3][i][j] = np.mean(cluster_Members)
                
        pickle.dump(clusterResults, open('sigScan.pickle','wb'))
        
    if (sys.argv[1] == '18'):
        """Plot Sig SCAN"""
        NMINRANGE = range(3,20,1) 
        signal = np.arange(.0,.4,.025)
        
        scale,S,count, members = pickle.load(open('epsScan.pickle','rb'))
        
        plt.figure(1)
        plt.subplot(221)
        CS = plt.contourf(S,10,extent = (signal[0],signal[-1],NMINRANGE[0], NMINRANGE[-1]))
        cbar = plt.colorbar(CS)
        cbar.ax.set_ylabel('S')
        plt.ylabel(r'$N_{min}$')
        plt.xlabel(r'$N_{sig}/N_{total}$')
        
        plt.subplot(222)
        CS = plt.contourf(scale,10,extent = (signal[0],signal[-1],NMINRANGE[0], NMINRANGE[-1]))
        cbar = plt.colorbar(CS)
        cbar.ax.set_ylabel('Scale')
        plt.ylabel(r'$N_{min}$')
        plt.xlabel(r'$N_{sig}/N_{total}$')
        
        plt.subplot(223)
        CS = plt.contourf(count,10,extent = (signal[0],signal[-1],NMINRANGE[0], NMINRANGE[-1]))
        cbar = plt.colorbar(CS)
        plt.ylabel(r'$N_{min}$')
        plt.xlabel(r'$N_{sig}/N_{total}$')
        cbar.ax.set_ylabel('Cluster Count')
        
        plt.subplot(224)
        CS = plt.contourf(members,10,extent = (signal[0],signal[-1],NMINRANGE[0], NMINRANGE[-1]))
        cbar = plt.colorbar(CS)
        cbar.ax.set_ylabel('Cluster Members')
        plt.ylabel(r'$N_{min}$')
        plt.xlabel(r'$N_{sig}/N_{total}$')
        
        plt.show()    
        
    #===========================================================================
    # SCAN EPS vs N_MIN
    #===========================================================================
    if (sys.argv[1] == '15'):
        profile = ('PULSAR',)
        fileOut = 'PULSARRateMap.pickle'
        
        # FERMI-LAT RUNS
        numPhotons2 = 48
        angularSize = 10.
        numSims = 200
        HESS = False
        for numPuls in [2,4,6]:
            MC.RUN_PULSAR(numSims, fileOut, numPhotons=numPhotons2,numPulsars=numPuls,angularSize=angularSize, outputSize=outputSize, mcList='test_' + str(numPuls) + '_' + str(numPhotons2)+ '.pickle',HESS=HESS,numProcs = numProcs)
        
        #HESS RUNS
        fileOut = 'PULSARRateMapHESS.pickle'
        HESS = True
        angularSize = 4.0
        numPhotons2 = 2000
        numSims = 10
        
        for numPuls in [2,4,6]:
            MC.RUN_PULSAR(numSims, fileOut, numPhotons=numPhotons2,numPulsars=numPuls,angularSize=angularSize, outputSize=outputSize, mcList='test_' + str(numPuls) + '_' + str(numPhotons2)+ '.pickle',HESS=HESS,numProcs = numProcs)
        

    
    if (sys.argv[1] == '13'):
        
        numProcs = 10
        angularSize = 10
        SNR = .75
        numAnalyze = 200
    
        HESS = False
        if HESS == True or sys.argv[2] == '2000':
            angularSize = 4.
            SNR = .95
            numAnalyze = 10
        
        mcSims = pickle.load(open('test_' + sys.argv[2] +'_' + sys.argv[3] + '.pickle', 'rb' ))
        
        NMINRANGE = range(3,10,1) 
        EPSRANGE = np.arange(.05,.40,.02)
        #NMINRANGE = range(4,8,1) 
        #EPSRANGE = np.arange(.08,.18,.02)
        size= len(NMINRANGE),len(EPSRANGE)
        clusterResults = [np.zeros(size),np.zeros(size),np.zeros(size),np.zeros(size)]
    
        for i in range(len(NMINRANGE)):
            for j in range(len(EPSRANGE)):
                #print 'eps, nmin: ', EPSRANGE[j], NMINRANGE[i]
                #cluster_Scale, cluster_S, cluster_Count, cluster_Members, stde v = [],[],[],[],[]
                cluster_Scale, cluster_S, cluster_Count, cluster_Members, stdev = MCSTATS.DBSCAN_Compute_Clusters(mcSims, S_cut=2.0, eps=EPSRANGE[j] , numAnalyze = numAnalyze,n_cluster=NMINRANGE[i], nCore = 3,HESS=HESS,angularSize = angularSize,SNR = SNR,numProcs=numProcs)
                #cluster_Scale, cluster_S, cluster_Count, cluster_Members = (DBSCAN_FERMI( energyRange = (120000,140000),eps = EPSRANGE[j], n_cluster = NMINRANGE[i], nCore = NMINRANGE[i]))
                clusterResults[0][i][j] = np.mean(cluster_Scale)
                clusterResults[1][i][j] = np.mean(cluster_S)
                clusterResults[2][i][j] = np.mean(cluster_Count)
                clusterResults[3][i][j] = np.mean(cluster_Members)
                print 'Progress: ', len(EPSRANGE)*i + j , '/' , len(NMINRANGE)*len(EPSRANGE)
                
        pickle.dump(clusterResults, open('epsScan_'+ sys.argv[2] +'_' + sys.argv[3] + '.pickle','wb'))
    

    if (sys.argv[1] == '14'):
        NMINRANGE = range(3,10,1) 
        EPSRANGE = np.arange(.05,.40,.02)
        
        #filenames = ['epsScan_2_48.pickle','epsScan_4_48.pickle','epsScan_6_48.pickle','epsScan_2_4000.pickle','epsScan_4_4000.pickle','epsScan_6_4000.pickle']
        
        fig2 = plt.figure(2, figsize = [12,8])
        
        epsPlot(1,'epsScan_2_48.pickle','2 Pulsar, Fermi-LAT')
        epsPlot(2,'epsScan_4_48.pickle','4 Pulsar, Fermi-LAT')
        epsPlot(3,'epsScan_6_48.pickle','6 Pulsar, Fermi-LAT')
        epsPlot(4,'epsScan_2_2000.pickle','2 Pulsar, HESS II')
        epsPlot(5,'epsScan_4_2000.pickle','4 Pulsar, HESS II')
        epsPlot(6,'epsScan_6_2000.pickle','6 Pulsar, HESS II')
        
        
        
        
#        plt.subplot(222)
#        CS = plt.contourf(scale,10,extent = (EPSRANGE[0],EPSRANGE[-1],NMINRANGE[0], NMINRANGE[-1]))
#        cbar = plt.colorbar(CS)
#        cbar.ax.set_ylabel('Scale')
#        plt.ylabel(r'$N_{min}$')
#        plt.xlabel(r'$\epsilon$')
#        
#        plt.subplot(223)
#        CS = plt.contourf(count,10,extent = (EPSRANGE[0],EPSRANGE[-1],NMINRANGE[0], NMINRANGE[-1]))
#        cbar = plt.colorbar(CS)
#        plt.ylabel(r'$N_{min}$')
#        plt.xlabel(r'$\epsilon$')
#        cbar.ax.set_ylabel('Cluster Count')
#        
#        plt.subplot(224)
#        CS = plt.contourf(members,10,extent = (EPSRANGE[0],EPSRANGE[-1],NMINRANGE[0], NMINRANGE[-1]))
#        cbar = plt.colorbar(CS)
#        cbar.ax.set_ylabel('Cluster Members')
#        plt.ylabel(r'$N_{min}$')
#        plt.xlabel(r'$\epsilon$')
        
        plt.savefig('eps_vs_nmin.pdf')
        plt.show()
        
    
    
    
    

#==============================================================================
# Radial Density Bins
#==============================================================================

from matplotlib import pyplot as plt
#mcSims1 = pickle.load(open('MCOut_galprop_NFW_500.pickle', "r" ))[0:2000]
#mcSims2 = pickle.load(open('MCOut_galprop_NFWDECAY_500.pickle', "r" ))[0:2000]

#NFW Rate Map
profile = ('NFW',23.5,1.0)
fileOut = 'NFWRateMap.pickle'
#mcSims1 = MC.RUN(100, fileOut, numPhotons=5000,angularSize=angularSize, outputSize=outputSize, mcList='MCOut_galprop_NFW_5000.pickle',HESS=True)
profile = ('NFWDECAY',23.5,1.0)
fileOut = 'NFWDECAYRateMap.pickle'
#mcSims2 = MC.RUN(100, fileOut, numPhotons=5000,angularSize=angularSize, outputSize=outputSize, mcList='MCOut_galprop_NFWDECAY_5000.pickle',HESS=True)
#mcSims3 = MC.RUN_BG_ONLY(100, fileOut, numPhotons=5000,angularSize=angularSize, outputSize=outputSize, mcList='MCOut_galprop_BG_5000.pickle',HESS=True)
#mcSims3 = MC.RUN_BG_ONLY(1000, fileOut, numPhotons=500,angularSize=angularSize, outputSize=outputSize, mcList='MCOut_galprop_BG_500.pickle',HESS=True)

#mcSims1 = pickle.load(open('MCOut_galprop_NFW_5000.pickle','rb'))
#mcSims2 = pickle.load(open('MCOut_galprop_NFWDECAY_5000.pickle','rb'))
#mcSims3 = pickle.load(open('MCOut_galprop_BG_5000.pickle','rb'))
#
#MCSTATS.Compute_Radial_Profiles(mcSims1,label = 'NFW 5e3')
#MCSTATS.Compute_Radial_Profiles(mcSims2,label = 'NFW Decay 5e3')
#MCSTATS.Compute_Radial_Profiles(mcSims3,label = 'Galprop BG 5e3')
#
#mcSims1 = pickle.load(open('MCOut_galprop_NFW_500.pickle','rb'))
#mcSims2 = pickle.load(open('MCOut_galprop_NFWDECAY_500.pickle','rb'))
#mcSims3 = pickle.load(open('MCOut_galprop_BG_500.pickle','rb'))
#
#MCSTATS.Compute_Radial_Profiles(mcSims1,label = 'NFW 500')
#MCSTATS.Compute_Radial_Profiles(mcSims2,label = 'NFW Decay 500')
#MCSTATS.Compute_Radial_Profiles(mcSims3,label = 'Galprop BG 500')
#
#plt.xlabel(r'$r[^\circ]$')
#plt.ylabel(r'$\rho [events deg^{-2}]$')
#plt.xlim(10**-1,10**1)
#plt.legend()
#plt.show()


