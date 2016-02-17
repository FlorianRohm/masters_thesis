from numpy import *

from centralMomentsFromDistributions import centralMomentsFromDistributions
from normalizedCumulantsFromCentralMoments import normalizedCumulantsFromCentralMoments
from normalizedCumulantsFromCentralMoments import normalizedCumulantsFromCentralMoments
from distributionsFromCentralMoments import distributionsFromCentralMoments

def normalizedCumulantsFromDistributions(rho, u, fin):
    (c_00, c_10, c_01, c_11, c_20, c_02, c_21, c_12, c_22) = centralMomentsFromDistributions (u, fin)

    return normalizedCumulantsFromCentralMoments (rho, u, c_00, c_10, c_01, c_11, c_20, c_02, c_21, c_12, c_22)

def distributionsFromNormalizedCumulants(rho, u, K_00, K_10, K_01, K_11, K_20, K_02, K_21, K_12, K_22):
    (c_00, c_10, c_01, c_11, c_20, c_02, c_21, c_12, c_22) = \
    normalizedCumulantsFromCentralMoments (rho, u, K_00, K_10, K_01, K_11, K_20, K_02, K_21, K_12, K_22)

    return distributionsFromCentralMoments (u, c_00, c_10, c_01, c_11, c_20, c_02, c_21, c_12, c_22)
