import numpy as np
from scipy.stats import kstest,uniform

def getperiodsampling(values, minsearchperiod, maxsearchperiod, oversample_factor=10):
    """Determine reasonable periods to sample when searching for mean period spacing.

    values are data to be assessed for mean spacing
    minsearchperiod is the shortest spacing to be considered
    maxsearchperiod is the longest spacing to be considered
    oversample_factor is how much more finely than critical sampling to sample periods
    """
    frequencyresolution = 1/(np.max(values)-np.min(values))
    frequencysample = np.arange(1/minsearchperiod,
                                1/maxsearchperiod-0.1*frequencyresolution/oversample_factor,
                                -frequencyresolution/oversample_factor)
    return 1/frequencysample #period


def wrapperiods(values, wrapperiod, centered=True):
    """Wrap values between 0-1 on a test period.

    values (array) to be wrapped on the number line from 0-1
    wrapperiod (float) to wrap the values on
    centered (bool; default:True) whether to center these around the "polar mean"
    """
    # Wrap on test period from 0-1
    mod = (values % wrapperiod)/wrapperiod
    
    if centered: #Center centroid on phase 0.5
        # Convert to radians to represent as vectors
        radians = mod*2*np.pi
        
        # Get x and y components of vectors
        xs = np.cos(radians)
        ys = np.sin(radians)
        
        # Get the average position vector components and angle
        meanx = np.average(xs) 
        meany = np.average(ys)
        meanangle = np.arctan2(meany,meanx)
        
        # How far is current average from being centered on pi radians?
        thetaoffset = (np.pi-meanangle) 
        
        # Shift everything by this offset
        radians  = (radians + thetaoffset) % (2*np.pi) # must put 2*np.pi in parenteses!
        
        # convert back to 0-1
        mod = radians/(2*np.pi)
    return mod


def IV(values, testperiods, centered = True):
    """Inverse Variance test

    values are data to be assessed for mean spacing
    test periods are spacings to evaluate
    centered says whether to center values on phase 0.5 (default: True, 
        you probably should)
    
    Returns the inverse variance of values for each test period
    """
    iv = np.zeros(len(testperiods))
    
    for i in range(len(testperiods)): #test different periods
        #mod test period
        wrapped = wrapperiods(values, testperiods[i], centered=centered)
        iv[i] = 1/np.var(wrapped)
    return iv
    

def KS(values, testperiods, centered = True):
    """Kolmogorov-Smirnov Test
    
    values are data to be assessed for mean spacing
    test periods are spacings to evaluate
    centered says whether to center values on phase 0.5 (default: True, 
        you probably should)
    
    Returns the log of the p-value (log Q) for each test period
    """
    Q = np.zeros(len(testperiods))
    
    for i in range(len(testperiods)): #test different periods
        #mod test period
        wrapped = wrapperiods(values, testperiods[i], centered=centered)
        Q[i] = kstest(wrapped,uniform.cdf).pvalue
    return np.log10(Q)


def FT(values, testperiods):
    """Fourier Transform Test
    
    values are data to be assessed for mean spacing
    test periods are spacings to evaluate
    
    Returns the power (best-fit zero-centerd sine amplitude squared) at each period tested
    """
    # Represent where there are period measurements as a window function
    window = np.ones(len(values))/2.
    
    # DFT function (stolen from Mikemon) fits zero-centered sine wave to values on test period
    ampvec = np.zeros(len(testperiods))
    frequencies = 1/testperiods
    for i, freq in enumerate(frequencies):
        omega = 2.*np.pi*freq
        wts = np.sin(omega*values)
        wtc = np.cos(omega*values)
        camp = np.dot(wtc, window)
        samp = np.dot(wts, window)
        ampvec[i] = np.sqrt(camp**2 + samp**2)
    ampvec = (2./len(values))*np.array(ampvec)
    
    return ampvec**2.
    