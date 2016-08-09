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

from numpy import roll
from getC import getC

c = getC()

def stream(f):
    for i in range(3):
        for j in range(3):
            for k in range(3):
                f[i,j,k] = roll(roll(roll( f[i,j,k], c[i,j,k,0], axis=0), c[i,j,k,1], axis=1), c[i,j,k,2], axis=2)
    return f
