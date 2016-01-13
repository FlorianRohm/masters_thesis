from centralMomentsFromCumulants import centralMomentsFromCumulants
from centralMomentsFromDistribution import centralMomentFromDistribution
from cumulantsFromCentralMoments import cumulantsFromCentralMoments
from distributionFromCentralMoments import distributionFromCentralMoments

def cumulantsFromDistribution(rhoInv, u, f):
    centralMoments = centralMomentFromDistribution(u, f);

    return cumulantsFromCentralMoments(rhoInv, *centralMoments)

def distributionFromCumulants(rhoInv, u, *cumulants):
    centralMoments = centralMomentsFromCumulants(rhoInv, *cumulants)

    return distributionFromCentralMoments(u, *centralMoments)
