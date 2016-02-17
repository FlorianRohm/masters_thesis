from numpy import *

from LBMHelpers import getMacroValues
from transformations.centralMomentsFromDistributions import centralMomentsFromDistributions
from transformations.distributionsFromCentralMoments import distributionsFromCentralMoments

def checkTransformation(f):
    (rho, u) = getMacroValues(f)

    centralMoments = centralMomentsFromDistributions(u, f)

    fNew = distributionsFromCentralMoments (u, *centralMoments)

    Delta_f = f-fNew

    print amax(Delta_f)
