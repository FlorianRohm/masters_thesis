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
from getWeights import getWeights



t = getWeights()

def equilibrium(rho, u):
    (d, nx, ny, nz) = u.shape
    cu = zeros((3, 3, 3, nx, ny, nz))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                cu[i, j, k] = (i-1)*u[0] + (j-1)*u[1] + (k-1)*u[2]
    usqr = 3./2.*(u[0]**2 + u[1]**2 + u[2]**2)
    feq = zeros((3, 3, 3, nx, ny, nz))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                feq[i, j, k] = rho*t[i,j,k]*(1. + 3.*cu[i,j,k] + 4.5*cu[i,j,k]**2 - usqr)
    return feq
