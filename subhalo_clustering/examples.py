#===============================================================================
# Examples of usage of Fermi Monte-Carlo Simulation Codes
# Author: Eric Carlson
# Updated: 11-14-2012
#===============================================================================

#===============================================================================
# Required Packages
#
#===============================================================================
# First of all, there are several packages that must be installed.
# These can be individually installed, or just download the 'Enthought Python 
# Distribution' which includes all of them.
#   You must have:  
#    numpy: http://numpy.scipy.org/
#    scipy: http://scipy.org/
#    sklearn: http://scikit-learn.org/stable/
#
# If you want to play around with the PSF, or build background templates from
# Fermi's galprop model you need pyfits to be included in the run folder.
#===============================================================================

#===============================================================================
# Module Descriptions:  Each module has comments and docstrings that describe
#     the parameters and methods.
#===============================================================================
# -BGTemplate.py: Generate background templates.  These are already included so you
#    do not need to do anything with this unless you want to play with the
#    background model.
# -Correlation.py: work in progress, but will compute the 2-pt correlation functions
#    using fermi data.
# -DBSCAN.py: Sort of a background module that runs DBSCAN and computes some of the
#     clustering properties.  For the most part you should just use MCSTATS.py. 
#    DBScan algorithm is in the sklearn package.
# -FermiPSF.py: Just converts the PSF fits from fermi tools into python tuples.
# -MC.py: This contains all simulation methods.
# -MCSTATS.py: Contains plotting functions, and wraps functionality of DBSCAN.py 
# -ParseFermi.py: Loads Fermi data from a csv file with specifications given in
#    the docstring. (photons.txt included has photons 10-300 GeV within 25 deg of GC)
# -runMC_v2.py has most of the code I have been using but it is in dissarray so just
#    use for reference if needed.
#===============================================================================

#===============================================================================
# Example usages
#===============================================================================
import pickle # serialization class
import sys, math

# Import the MC class that has our simulation code.
import MC, MCSTATS

#First let us specify some global parameters.
angularSize = 10.0 # width of box centered on galactic center
outputSize = 300  # grid size for background template, ratemap, etc.  Not really important as long as angularSize/outputSize> instrument resolution.
numSims = 100 # Number of simulations to run.
numPhotons = 48 # number of photons in a simulation

# Now let's generate the J-factor (LOS integral of annihilation) for an NFW profile of scaling radius 23.5 and alpha=1
# Note that the profile specifications can be found in the MC.RUN docstrings 
profile = ('NFW',23.5,1.0) 



# Specify where to store the rate map.  This only needs to be generated once for a given profile and then we can just reference the file.  
# I have included the distributions we have already done in dropbox  
fileOut = 'NFWRateMap.pickle' # location to save the output projected annihilation map
#MC.Gen_Annihilation_Map(angularSize, outputSize, profile,fileOut)

# If we want to view it we can uncomment the following lines and it will save a png too named fileOut.png
#MCSTATS.Plot_Rate_Map('NFWRateMap.pickle', angularSize,fileOut)

# Now let's run a set of simulations.  The HESS parameter forces a tighter PSF to be used instead of distributing between back and front converting events
mcSims = MC.RUN(numSims,rateMap=fileOut, numPhotons=numPhotons,angularSize=angularSize, outputSize=outputSize, mcList='MCOut_galprop_NFW_48_example.pickle',HESS=False)

# Now our simulations are stored in the simulation list 'mcSims' object, or the pickle file so we don't have to run simulations again. 
# If we like we can plot the first simulation using:
#MCSTATS.Plot_MC_Positions(MC=mcSims[0], angularSize=angularSize)

# Let's compute a pulsar simulation set to so we can compare.  This is a slightly different method, but almost identical
profile = ('PULSAR',)
fileOut = 'PULSARRateMap.pickle'
# Generate rate map if not done already, but should be included.
#MC.Gen_Annihilation_Map(angularSize, outputSize, profile,fileOut)

# Run the pulsar version of the MC
numPulsars = 3
mcSims = MC.RUN_PULSAR(numSims, fileOut, numPhotons=numPhotons,numPulsars=numPulsars,angularSize=angularSize, outputSize=outputSize, mcList='MCOut_galprop_PULSAR_48_' + str(numPulsars) + '_example.pickle')



# Presumably we want to get some statistics out so let's load the simulation results from file and run DBSCAN 
# (most of these parameters are default values, but shown here explicitly)

NFWSims = pickle.load(open('MCOut_galprop_NFW_48_example.pickle','rb'))
pulsarSims = pickle.load(open('MCOut_galprop_PULSAR_48_' + str(numPulsars) + '_example.pickle','rb'))

eps=0.1875         # 95% PSF containment for rear-converting events
MCSTATS.DBSCAN_Compute_Clusters(mcSims=NFWSims, eps=eps, n_cluster=3, nCore=3, S_cut=1.0, numAnalyze=0, fileout='Clusters_NFW.pickle')
MCSTATS.DBSCAN_Compute_Clusters(mcSims=pulsarSims, eps=eps, n_cluster=3, nCore=3, S_cut=1.0, numAnalyze=0, fileout='Clusters_pulsar.pickle')
# Now our clusters are stored in Clusters_*.pickle so we will not have to run clustering again on this simulation set.

# Finally, we want to compute statistics. so let us define a function for all the matplotlib details.
import matplotlib.pyplot as plt 
def Gen_Plots(fileout = 'Cluster_Example_48'):
    """Compare Simulations"""
    inputList = (('NFW','Clusters_NFW.pickle'),
                 ('Pulsar 3','Clusters_pulsar.pickle'))
                 
    labels = []
    scales = []
    sigs   = []
    count  = []
    members= []
    for i in inputList:
        (cluster_Scale, cluster_S, cluster_Count, cluster_Members,cluster_stdevs) = pickle.load(open(i[1],'r'))
        scales.append(cluster_Scale)
        sigs.append(cluster_S)
        count.append(cluster_Count)
        members.append(cluster_Members)
        labels.append(i[0])
    
    fig = plt.figure(1, (8,14))
    plt.clf()
    
    from matplotlib.font_manager import FontProperties
    fontP = FontProperties()
    fontP.set_size('small')
    
    MCSTATS.Plot_Cluster_Scales(sigs, labels,r'S',fig,1,   bins=50, fileout='Cluster_Sigs',PlotFermiScale = True)
    plt.legend(loc=1, ncol=1, fancybox=True, shadow=False,prop=fontP,borderaxespad=0.,labelspacing = .2)
    MCSTATS.Plot_Cluster_Scales(scales, labels,r'$r_{cluster}[^\circ]$',fig,2, bins=50, fileout='Cluster_Scaling',PlotFermiScale = True)
    MCSTATS.Plot_Cluster_Scales(count, labels,r'$N_{clusters} | S>2.0$',fig,3,  bins=50, fileout='Cluster_Count',PlotFermiScale = True)
    MCSTATS.Plot_Cluster_Scales(members, labels,r'$N_{members}$',fig,4,bins=50, fileout='Cluster_Members',PlotFermiScale = True)
    
    if fileout != '':
        from matplotlib.backends.backend_pdf import PdfPages
        pp = PdfPages(fileout + '.pdf')
        plt.savefig(pp, format='pdf')
        print "Figures saved to ", str(fileout)+ '.pdf\n',
        pp.close()

    return
 
# Call function and show plot.
Gen_Plots()
plt.show()



#===============================================================================
# That is about it.  The rest can be inferred from runMC_v2, which also has some
# details of how the fermi data is accessed in runMC_v2.DBSCAN_FERMI()
#===============================================================================






