from numpy import *

from centralMomentsFromDistributions import centralMomentsFromDistributions
from normalizedCumulantsFromCentralMoments import normalizedCumulantsFromCentralMoments
from centralMomentsFromNormalizedCumulants import centralMomentsFromNormalizedCumulants
from distributionsFromCentralMoments import distributionsFromCentralMoments

from centralMomentsFromDistributions_min import centralMomentsFromDistributions_min
from centralMomentsFromNormalizedCumulants_min import centralMomentsFromNormalizedCumulants_min
from distributionsFromCentralMoments_min import distributionsFromCentralMoments_min

def normalizedCumulantsFromDistributions(rho, u, fin):
    centralMoments = centralMomentsFromDistributions (u, fin)

    return normalizedCumulantsFromCentralMoments (rho, u, *centralMoments)

def distributionsFromNormalizedCumulants(rho, u, *normalizedCumulants):
    centralMoments = centralMomentsFromNormalizedCumulants (rho, *normalizedCumulants)

    return distributionsFromCentralMoments (u, *centralMoments)

def normalizedCumulantsFromDistributions_min(rho, u, fin):
    return centralMomentsFromDistributions_min (u, fin)

def distributionsFromNormalizedCumulants_min(rho, u, *normalizedCumulants):
    centralMoments = centralMomentsFromNormalizedCumulants_min (rho, *normalizedCumulants)

    return distributionsFromCentralMoments_min (u, *centralMoments)
