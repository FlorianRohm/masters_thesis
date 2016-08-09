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

from numpy import zeros

def getC():
    c = zeros((3, 3, 3, 3), dtype='int64')
    # dimensions of c are
    # the three indices of D3Q27, ijk
    # the velocity vector for the index
    for i in range(3):
        for j in range(3):
            for k in range(3):
                c[i, j, k, 0] = (i-1)
                c[i, j, k, 1] = (j-1)
                c[i, j, k, 2] = (k-1)
    return c
