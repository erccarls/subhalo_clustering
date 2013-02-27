#===============================================================================
# MC.py:  Contains all Monte-Carlo methods for simulating various source classes
#  in presence of a continuum backgound.
# Author: Eric Carlson 
# Updated: 11-14-2012
#===============================================================================
import numpy as np
import math,pickle,sys,time  

import multiprocessing as mp
from multiprocessing import pool
from functools import partial


def Gen_Annihilation_Map(angularSize, size, profile,fileOut):
    """
    Generate map of the annihilation rate (i.e. LOS integral), centered on the galactic 
    center based on an input DM profile. The rate is scaled such that max = 1.  In the 
    case of a pulsar profile this will build a map of the probability to place a pulsar
    at a given lat. lon., but the actual puslar placements are randomly chosen 
    accordingly at a later time.  Generally this should just be generated once and then 
    saved for all future use, but it does return the map directly as well. 
    
    Input: 
        angularSize: Number of degrees to span
        size: Number of pixels per side (square)
        profile: ('NFW', rc, alpha), ('EIN',rc,alpha),('NFWDECAY', rc, alpha), ('FLAT',0,0),('PULSAR',0,0):
        fileOut: Output pickle filename for rate map
        
    Return: 
        -rateMap: Array of relative annihilation rates normalized to max = 1.
    """
    print 'Generating Annihilation Rate Map'
    
    # Behavioral Parameters    
    solarDist = 8.3                        # Solar Distance in kpc
    stopRadius = 60.0                      # Radius from sun in kpc to stop LOS integration
    zSteps = 100                           # Steps for LOS integration
    kpc2cm = 3.08568025e21
    
    # Constants
    deg2rad = math.pi/180.0 
    map = np.zeros((size, size))             # Initialize map
    
    # Code
    zStepSize = stopRadius/float(zSteps)                 # Size of one step in kpc
    aPP = float(angularSize)/float(size)                 # Angle per pixel
    solidAngle = (aPP*deg2rad)**2.0                      # Multiply this by radial portion to get wedge volume.
    # Based on the input DM_model, we integrate rho**2 along the LOS
    max = 0
    for x in range(0,size):               
        for y in range(0,size):
            rate = 0.0
            gamma = math.sqrt((float(x)-size/2.0)**2 + (float(y)-size/2.0)**2)*aPP # Inclusive angle for law of cosines
            if profile[0] != 'PULSAR': 
                for z in range(0,zSteps):
                    # Compute wedge volume.  Currently this assumes a relatively small angular region around galactic center.
                    volume = ((z+1)**3.0-z**3.0)*solidAngle/3.#
                    # Compute radius from galactic center using law of cosines
                    a = (z*zStepSize)
                    #r = math.sqrt(a**2 + b**2 - 2*a*b*math.cos(gamma)) 
                    l = (float(x)-size/2.0)*aPP # longitude
                    b = (float(y)-size/2.0)*aPP # latitude 
                    r = math.sqrt(a**2 + solarDist**2 - 2*a*solarDist*math.cos(l*deg2rad)*math.cos(b*deg2rad))
                    #if gamma>max:
                    #    max = gamma
                    # Get square DM density to obtain rate
                    if (profile[0] == 'NFW'):
                        rate += (volume*rho_DM_NFW(r,profile[1],profile[2]))**2
                    elif (profile[0] == 'EIN'):
                        rate += (0.0780763*rho_DM_EIN(r,profile[1],profile[2]))**2*zStepSize*kpc2cm*solidAngle
                        
                    elif (profile[0] == 'FLAT'):
                        rate = 1 # just keep everything flat.  += will give integration limit dependent results
                    elif (profile[0] == 'NFWDECAY'): # NFW not squared.
                        rate += volume*rho_DM_NFW(r,profile[1],profile[2])
                map[x,y] = rate
            else:
                l = (float(x)-size/2.0)*aPP # longitude
                b = (float(y)-size/2.0)*aPP # latitude
                r = math.sqrt(l**2.0+b**2.0)
                if (r<=0.05):
                    map[x,y] = 0.05**-1.6
                else:
                    map[x,y] = r**-1.2 
    # Write to file        
    outFile = open(fileOut, "wb" )
    pickle.dump(map/np.max(map), outFile)
    print 'Rate Map saved to ', fileOut
    #print np.max(map)
    print 'J-Factor (GeV^2/cm^5): ' , np.sum(map)
    return map/np.max(map)


def rho_DM_NFW(r,rc,alpha):
    if (r != 0.0):
        return (rc/r)**alpha * (1+ r/rc)**(-3+alpha)
    else:
        return 0.0

def rho_DM_EIN(r,rc,alpha):
    return math.exp(-2.0/alpha*((r/rc)**alpha-1.))


def RUN(numTrials, rateMap, numPhotons=48, angularSize=10.0, outputSize=300, mcList='MCOut.pickle',HESS=False, Sig = -1 ,numProcs = 10):
    """
    Runs a set of MC simulations for diffuse emission sources and outputs a list of results for each simulation 
    (IN PIXEL COORDINATES: Must be shifted and scaled from pixels -> degrees if using these results directly. This 
    is already done in all written methods)
    
    Inputs:
        -numTrials: number of simulations to run
        -rateMap: the pickled array output from Gen_Annihilation_Map()
        -numPhotons: number of photons in each simulation.
        -angularSize: Size of simulation in degrees
        -outputSize: Size of simulation in pixels.  Needs to match rateMap.
        -mcList: output results filename
        -HESS: true uses only fermi's front converting PSF since this more closely matches ACT's.
               false distributes events conversion type randomly according to effective exposure
               area of front vs. back.
        -Sig, if a custom Signal percentage is wanted, pass it here.  otherwise 25% is used for Fermi and 5% for HESS
               
    Returns:
        MCOut: list or coordinate pairs of points in degrees.
    """
    print 'Beginning MC Series\nProgress'
    
    import FermiPSF, ParseFermi
    mcOut = []
    map = pickle.load(open(rateMap, "r" )) # load rate-map
    PSFTableFront = FermiPSF.PSF_130(convType='front') # load PSF front converting
    PSFTableBack = FermiPSF.PSF_130(convType='back')   # load PSF back converting

    start = time.time();
    
    ppa = outputSize/angularSize # pixel per degree

    # Import background template
    bgmap = 'BGRateMap.pickle'
    if (HESS == True):
        bgmap = 'BGRateMap_HESS_2_deg.pickle'
    
    bgTemplate = pickle.load(open(bgmap , "r" ))
    
    mcOut = np.zeros(numTrials)
    p = pool.Pool(numProcs)
    
    partial_MC_THREAD = partial( MC_THREAD, map = map,bgTemplate=bgTemplate,PSFTableFront=PSFTableFront, PSFTableBack=PSFTableBack, HESS=HESS, angularSize=angularSize, numPhotons=numPhotons, outputSize=outputSize,Sig = Sig)
    mcOut = p.map(partial_MC_THREAD, mcOut)
    
#    for i in range(numTrials):   
#        # Build the background                    
##        background = Build_Background_Sideband(bgMean, lowSideband, highSideband, PSFTable)
#        background = Build_Background_Template(bg, bgTemplate, PSFTableFront, PSFTableBack,flatLevel = 0.0,HESS= HESS,angularSize = angularSize)
#        # Compute number of source photons
#        numMC  = numPhotons - len(background[0])
#        # Run MC for source photons       
#        data = MC(map,numMC,angularSize,outputSize,PSFTableFront, PSFTableBack,HESS=HESS)
#        # Append data
#        mcOut.append((data[0]+background[0], data[1]+background[1]))
#        
#        # Compute Speed Statistics
#        sys.stdout.write('\r' + str(i+1)+'/'+str(numTrials)) 
#        sys.stdout.flush()
    elapsed = time.time()-start;
    if (elapsed != 0.0):
        print '\nSimulations Completed in', elapsed, 's', '(',numTrials/elapsed, ' sims per second)'
    
    outFile = open(mcList, "wb" )
    pickle.dump(mcOut, outFile)
    print 'Results saved to ', mcList
    return mcOut


def MC_THREAD(x, map,bgTemplate,PSFTableFront, PSFTableBack, HESS, angularSize, numPhotons, outputSize,Sig):
       
        np.random.seed()
        # Compute number of background photons
        numSignal = np.random.poisson(lam = .25*numPhotons)
        if (HESS == True):
            numSignal = np.random.poisson(lam = .05*numPhotons)
        if Sig >= 0:
            numSignal = np.random.poisson(lam = Sig*numPhotons)                
        # Build the background 
        bg = numPhotons-numSignal # number of BG photons
        background = Build_Background_Template(bg, bgTemplate, PSFTableFront, PSFTableBack,flatLevel = 0.0,HESS= HESS,angularSize = angularSize)
        # Compute number of source photons
        numMC  = numPhotons - len(background[0])
        # Run MC for source photons       
        data = MC(map,numMC,angularSize,outputSize,PSFTableFront, PSFTableBack,HESS=HESS)
        # Append data
        return (data[0]+background[0], data[1]+background[1])

def MC_PULSAR_THREAD(x, map,bgTemplate,PSFTableFront, PSFTableBack, HESS, angularSize, numPhotons, outputSize,numPulsars,Sig):
    
        np.random.seed()
        # Compute number of background photons
        numSignal = np.random.poisson(lam = .25*numPhotons)
        if (HESS == True):
            numSignal = np.random.poisson(lam = .05*numPhotons)
        if Sig >= 0:
            numSignal = np.random.poisson(lam = Sig*numPhotons)
        
        # Build the background                    
        bg = numPhotons-numSignal # number of BG photons
        background = Build_Background_Template(bg, bgTemplate, PSFTableFront, PSFTableBack ,HESS=HESS, angularSize = angularSize )
        # Compute number of source photons
        numMC  = numPhotons - len(background[0])
        
        # Run MC for source photons       
        data = MC_PULSAR(map,numMC, numPulsars,angularSize,outputSize,PSFTableFront, PSFTableBack, HESS = HESS)
        # Concatenate and append this run to the simulation output
        return(data[0]+background[0], data[1]+background[1])


    
def RUN_PULSAR(numTrials, rateMap, numPhotons=48,numPulsars = 6, angularSize=10.0, outputSize=100, mcList='MCOut.pickle',flatLevel = 0.0,HESS=False, Sig = -1,numProcs = 10):
    """
    Runs numTrials MC simulations and outputs a pickle file with concatenated results ((x1,y1),(x2,y2),...) 
    """
    import FermiPSF, ParseFermi
    
    print 'Beginning MC Series\nProgress'

    mcOut = []
    map = pickle.load(open(rateMap, "r" )) # load rate-map
    PSFTableFront = FermiPSF.PSF_130(convType='front') # load PSF front converting
    PSFTableBack = FermiPSF.PSF_130(convType='back')   # load PSF back converting
    start = time.time();
     
    ppa = outputSize/angularSize # pixel per degree

    # Import background template
    bgmap = 'BGRateMap.pickle'
    if (HESS == True):
        bgmap = 'BGRateMap_HESS_2_deg.pickle'
        
    bgTemplate = pickle.load(open(bgmap , "r" ))
    
    mcOut = np.zeros(numTrials)
    p = pool.Pool(numProcs)
    partial_MC_PULSAR_THREAD = partial( MC_PULSAR_THREAD, map = map,bgTemplate=bgTemplate,PSFTableFront=PSFTableFront, PSFTableBack=PSFTableBack, HESS=HESS, angularSize=angularSize, numPhotons=numPhotons, outputSize=outputSize, numPulsars = numPulsars,Sig=Sig)
    mcOut = p.map(partial_MC_PULSAR_THREAD, mcOut)
    
#    for i in range(numTrials):
#        np.random.seed()
#        # Compute number of background photons
#        numSignal = np.random.poisson(lam = .25*numPhotons)
#        if (HESS == True):
#            numSignal = np.random.poisson(lam = .05*numPhotons)
#        if Sig >= 0:
#            numSignal = np.random.poisson(lam = Sig*numPhotons)
#        
#        bg = numPhotons-numSignal # number of BG photons
#   
#        # Build the background                    
##        background = Build_Background_Sideband(bgMean, lowSideband, highSideband, PSFTable)
#        background = Build_Background_Template(bg, bgTemplate, PSFTableFront, PSFTableBack ,HESS=HESS, angularSize = angularSize )
#        
#        
#        # Run MC for source photons       
#        data = MC_PULSAR(map,numSignal, numPulsars,angularSize,outputSize,PSFTableFront, PSFTableBack, HESS = HESS)
#        # Concatenate and append this run to the simulation output
#        mcOut.append((data[0]+background[0], data[1]+background[1]))
#        
#        # Compute Speed Statistics
#        sys.stdout.write('\r' + str(i+1)+'/'+str(numTrials)) 
#        sys.stdout.flush()
    elapsed = time.time()-start;
    if (elapsed != 0.0):
        print '\nSimulations Completed in', elapsed, 's', '(',numTrials/elapsed, ' sims per second)'
    
    outFile = open(mcList, "wb" )
    pickle.dump(mcOut, outFile)
    print 'Results saved to ', mcList
    return mcOut    


def RUN_BG_ONLY(numTrials, rateMap, numPhotons=48, angularSize=10.0, outputSize=300, mcList='MCOut.pickle',HESS=False):
    """
    Runs a set of MC simulations of just a background model  and outputs a list of results for each simulation 
    (IN PIXEL COORDINATES: Must be shifted and scaled from pixels -> degrees if using these results directly. This 
    is already done in all written methods)
    
    Inputs:
        -numTrials: number of simulations to run
        -rateMap: the pickled array output from Gen_Annihilation_Map() (UNUSED HERE)
        -numPhotons: number of photons in each simulation.
        -angularSize: Size of simulation in degrees
        -outputSize: Size of simulation in pixels.  Needs to match rateMap.
        -mcList: output results filename
        -HESS: true uses only fermi's front converting PSF since this more closely matches ACT's.
               false distributes events conversion type randomly according to effective exposure
               area of front vs. back.
               
    Returns:
        MCOut: list or coordinate pairs of points in degrees.
    """
    print 'Beginning MC Series\nProgress'
    
    import FermiPSF, ParseFermi
    mcOut = []
    map = pickle.load(open(rateMap, "r" )) # load rate-map
    PSFTableFront = FermiPSF.PSF_130(convType='front') # load PSF front converting
    PSFTableBack = FermiPSF.PSF_130(convType='back')   # load PSF back converting

    start = time.time();
    
    # Compute number of background photons
    bgMean = numPhotons # number of BG photons
    ppa = outputSize/angularSize # pixel per degree

    # Import background template
    bgTemplate = pickle.load(open('BGRateMap.pickle' , "r" ))
    
    for i in range(numTrials):   
        # Build the background                    
#        background = Build_Background_Sideband(bgMean, lowSideband, highSideband, PSFTable)
        background = Build_Background_Template(bgMean, bgTemplate, PSFTableFront, PSFTableBack,flatLevel = 0.0,HESS= HESS)
        # Compute number of source photons
        numMC  = numPhotons - len(background[0])
        # Append data
        mcOut.append((background[0], background[1]))
        
        # Compute Speed Statistics
        sys.stdout.write('\r' + str(i+1)+'/'+str(numTrials)) 
        sys.stdout.flush()
    elapsed = time.time()-start;
    if (elapsed != 0.0):
        print '\nSimulations Completed in', elapsed, 's', '(',numTrials/elapsed, ' sims per second)'
    
    outFile = open(mcList, "wb" )
    pickle.dump(mcOut, outFile)
    print 'Results saved to ', mcList
    return mcOut
    
def MC(rateMap,numPhotons,angularSize,outputSize,PSFTableFront, PSFTableBack, HESS=False):
    """
    Given a specified number of events, this uses a map of annihilation rate along the LOS to generate events.  A point spread function is then used for each event to determine the location of the simulated event in the instrument.
    Input:
        rateMap:    A square array scaled from 0-1 that specifies relative annihilation rate along LOS
        numPhotons: Number of photons to produce in the MC
        angularSize: AngularSize of the input rateMap
        outputSize:  Size of the output map.  This should correspond to the number of pixels from Fermi that span the angularSize
    Output: 
        Returns a list of x coordinates and of y coordinates for generated events.
    """
    
    x_dimen = np.shape(rateMap)[0]
    y_dimen = np.shape(rateMap)[1]
    
    outputScaleFactor = float(outputSize)/float(x_dimen)
    APP = float(angularSize)/float(outputSize)
    
    photonListX = []
    photonListY = []
    
    photonCount = 0
    while (photonCount <numPhotons):
        # Choose a random coordinate in the rate map
        x = np.random.randint(0,high=x_dimen)
        y = np.random.randint(0,high=y_dimen)
        # Look up value of annihilation rate
        rate = rateMap[x,y]
        # Select random number between 0 and 1.  If rate is greater, then we accept an annihilation here
        if (np.random.ranf() < rate):
            # Shift and scale coordinates to output map and then compute PSF modification to the position.
            psfMod = PSF_Spread(PSFTableFront,PSFTableBack, HESS = HESS)
            dx = psfMod[0]*math.cos(psfMod[1])/APP # PSF shift in output pixels
            dy = psfMod[0]*math.sin(psfMod[1])/APP # 
            x = int(round((x)*outputScaleFactor + dx)) 
            y = int(round((y)*outputScaleFactor + dy))
            # Ensure that we are still in the region of interest after PSF modification
            if (abs(x) <= outputSize and abs(y) <= outputSize):
                photonListX.append((x-float(outputSize)/2.0)*APP)
                photonListY.append((y-float(outputSize)/2.0)*APP)
                photonCount+=1  
    return (photonListX,photonListY)


def MC_PULSAR(rateMap,numPhotons,numPulsars,angularSize,outputSize,PSFTableFront, PSFTableBack, HESS = False):
    """
    Given a specified number of events, this uses a map of annihilation rate along the LOS to generate events.  A point spread function is then used for each event to determine the location of the simulated event in the instrument.
    Input:
        rateMap:    A square array scaled from 0-1 that specifies relative annihilation rate along LOS
        numPhotons: Number of photons to produce in the MC
        angularSize: AngularSize of the input rateMap and output region
        outputSize:  Size of the output map.  This should correspond to the number of pixels from Fermi that span the angularSize
    Output: 
        Returns a list of x coordinates and of y coordinates for generated events.
    """
    
    x_dimen = np.shape(rateMap)[0]
    y_dimen = np.shape(rateMap)[1]
    
    outputScaleFactor = float(outputSize)/float(x_dimen)
    
    APP = float(angularSize)/float(outputSize)
    
    photonListX = []
    photonListY = []
    photonCount = 0
    
    pulsarListX = []
    pulsarListY = []
    pulsarCount = 0
    
    # First pick out the projected position of the pulsars
    while (pulsarCount < numPulsars):
        # Choose a random coordinate in the rate map
        x = np.random.randint(0,high=x_dimen)
        y = np.random.randint(0,high=y_dimen)
        # Look up value of annihilation rate
        rate = rateMap[x,y]
        # Select random number between 0 and 1.  If rate is greater, then we accept a pulsar here
        if (np.random.ranf() < rate):
            # Shift and scale coordinates to output map
            x = int(round((x)*outputScaleFactor)) 
            y = int(round((y)*outputScaleFactor))
            # Ensure that we are still in the region of interest after PSF modification
            if (abs(x) <= outputSize and abs(y) <= outputSize):
                pulsarListX.append(x)
                pulsarListY.append(y)
                pulsarCount+=1
                
    # Now for each photon we must choose a progenitor pulsar and modify the photon position by the Fermi PSF
    while (photonCount <numPhotons):
        # Choose a random pulsar (need to weight by distance?)
        idx = np.random.randint(0,high=numPulsars)
        
        x = pulsarListX[idx]
        y = pulsarListY[idx]
        
        # Currently equal weight given to each pulsar so no need for anything but PSF modification.
    
        # Shift and scale coordinates to output map and then compute PSF modification to the position.
        psfMod = PSF_Spread(PSFTableFront, PSFTableBack, HESS = HESS)
        dx = psfMod[0]*math.cos(psfMod[1])/APP # PSF shift in output pixels
        dy = psfMod[0]*math.sin(psfMod[1])/APP # 
        x = x*outputScaleFactor + dx 
        y = y*outputScaleFactor + dy
        
        # Ensure that we are still in the region of interest after PSF modification
        if (abs(x) <= outputSize and abs(y) <= outputSize):
            photonListX.append((x-float(outputSize)/2.0)*APP)
            photonListY.append((y-float(outputSize)/2.0)*APP)
            photonCount+=1  
    return (photonListX,photonListY)


def Build_Background_Template(numBGPhotons, bgTemplate, PSFTableFront, PSFTableBack,flatLevel = 0.0,HESS = False,outputSize=300,angularSize=10.0):
    """
    Builds a background template based on sidebands.
        
        bgtemplate should be a array with a background model normalized to max(bgTemplate) = 1.0
        APP is the angle per pixel in degrees
    """
    
    numPhotons = numBGPhotons
    numHigh = int(round(.32 *numPhotons))
    numLow  = numPhotons-numHigh
    
    bgEventsX = []
    bgEventsY = []
    
    bgTemplate = bgTemplate *(1.0-flatLevel) + flatLevel
#    import matplotlib.pyplot as plt
#    plt.imshow(bgTemplate,'jet',vmin=0, vmax=1)
#    plt.colorbar()
#    plt.show()

    app=float(angularSize)/float(outputSize) # angle per pixel
    for i in range(numPhotons):
        x ,y = 0, 0
        while True:
            x,y = np.random.randint(0,high = len(bgTemplate)),np.random.randint(0,high = len(bgTemplate))
            if (np.random.ranf() < bgTemplate[y][x]):
                break
        # Shift and scale coordinates to output map and then compute PSF modification to the position.
        psfMod = PSF_Spread(PSFTableFront,PSFTableBack, HESS =HESS)
        dx = psfMod[0]*math.cos(psfMod[1]) # PSF shift in deg
        dy = psfMod[0]*math.sin(psfMod[1]) # PSF shift in deg
        
        bgEventsX.append((x-outputSize/2.0)*app + dx)
        bgEventsY.append((y-outputSize/2.0)*app + dy)
        
    return (bgEventsX, bgEventsY) 


def PSF_Spread(psffront,psfback,HESS = False):
    radialCutoff = 1.0    # Not going to be more than ... degrees
    r = 0.0
    psf = psffront
    
    FrontArea = 0.561 # percentage of front converting event.  The rest are back converting.
    if ((np.random.ranf() >= FrontArea) and (HESS == False)):
        psf = psfback
    
    psfMax = np.max(psf[1])
    while (True):
        r = np.random.ranf()*radialCutoff # scale this to cutoff radius.
        idx=(np.abs(psf[0]-r)).argmin() 
        distVal = psf[1][idx]
        y = np.random.ranf()*psfMax
        if (y<distVal):
            break    
    theta = np.random.ranf()*2*math.pi
    return (r,theta)







    
    