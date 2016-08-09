# 2D Lattice Boltzmann Code
# Copyright (C) 2016  Florian Rohm
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
