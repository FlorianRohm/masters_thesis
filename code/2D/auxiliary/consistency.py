from numpy import *

from LBMHelpers import getMacroValues
from transformations.centralMomentsFromDistributions import centralMomentsFromDistributions
from transformations.distributionsFromCentralMoments import distributionsFromCentralMoments
from transformations.normalizedCumulantsFromCentralMoments import normalizedCumulantsFromCentralMoments
from transformations.centralMomentsFromNormalizedCumulants import centralMomentsFromNormalizedCumulants
from transformations.transforms import normalizedCumulantsFromDistributions, distributionsFromNormalizedCumulants

def checkTransformation(f):
    threshold = 1./100000000.
    (rho, u) = getMacroValues(f)

    centralMoments = centralMomentsFromDistributions(u, f)
    fNew = distributionsFromCentralMoments (u, *centralMoments)

    Delta_f = f-fNew
    if amax(Delta_f) > threshold:
        print amax(Delta_f)

    centralMoments = centralMomentsFromDistributions(u, f)
    cumulants = normalizedCumulantsFromCentralMoments(rho, u, *centralMoments)
    centralMomentsNew = centralMomentsFromNormalizedCumulants(rho, *cumulants)
    fNew = distributionsFromCentralMoments (u, *centralMoments)

    Delta_f = f-fNew
    if amax(Delta_f) > threshold:
        print amax(Delta_f)

    cumulants = normalizedCumulantsFromDistributions(rho, u, f)
    fNew = distributionsFromNormalizedCumulants(rho, u, *cumulants)

    Delta_f = f-fNew
    if amax(Delta_f) > threshold:
        print amax(Delta_f)
