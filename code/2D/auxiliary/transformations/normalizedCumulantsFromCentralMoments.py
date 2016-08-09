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

def normalizedCumulantsFromCentralMoments (rho, u, c_00, c_10, c_01, c_11, c_20, c_02, c_21, c_12, c_22):

    ux = u[0, :, :]
    uy = u[1, :, :]

    K_22 = c_22 - 2*c_11*c_11/rho - c_20*c_02/rho

    return (log(rho)*rho, ux*rho, uy*rho, c_11, c_20, c_02, c_21, c_12, K_22)
