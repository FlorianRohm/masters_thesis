from numpy import *

from centralMomentsFromDistributions import centralMomentsFromDistributions
from normalizedCumulantsFromCentralMoments import normalizedCumulantsFromCentralMoments
from centralMomentsFromNormalizedCumulants import centralMomentsFromNormalizedCumulants
from distributionsFromCentralMoments import distributionsFromCentralMoments

def normalizedCumulantsFromDistributions(rho, u, fin):
    centralMoments = centralMomentsFromDistributions (u, fin)

    return normalizedCumulantsFromCentralMoments (rho, u, *centralMoments)

def distributionsFromNormalizedCumulants(rho, u, *normalizedCumulants):
    centralMoments = centralMomentsFromNormalizedCumulants (rho, *normalizedCumulants)

    return distributionsFromCentralMoments (u, *centralMoments)
