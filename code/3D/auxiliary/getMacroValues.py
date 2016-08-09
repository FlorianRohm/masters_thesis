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

from numpy import *
from getC import getC

c = getC()
def sumPopulations(fin):
    return sum(fin, axis = (0,1,2))

def getMacroValues(f):
    (d1,d2,d3, nx, ny, nz) = f.shape
    u = zeros((3, nx, ny, nz))

    # Calculate macroscopic density ...
    rho = sumPopulations(f)
    # ... and velocity
    zeros((3, 3, 3, nx, ny, nz))
    for i in range(3):
        for j in range(3):
            u[0] = u[0] + c[2, i, j, 0] * f[2, i, j] + c[0, i, j, 0] * f[0, i, j]
            u[1] = u[1] + c[i, 2, j, 1] * f[i, 2, j] + c[i, 0, j, 1] * f[i, 0, j]
            u[2] = u[2] + c[i, j, 2, 2] * f[i, j, 2] + c[i, j, 0, 2] * f[i, j, 0]
    u = u/rho
    return (rho, u)
