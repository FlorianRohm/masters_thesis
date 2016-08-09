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
