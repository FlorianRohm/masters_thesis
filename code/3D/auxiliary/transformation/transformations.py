# 3D Lattice Boltzmann Code
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
