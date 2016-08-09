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

from pyevtk.hl import gridToVTK

def saveToVTK(rho, u, prefix, saveNumber, grid):
    name = "./" + prefix + "." + saveNumber
    reorderedU = (u[0], u[1], u[2])
    gridToVTK(name, grid[0], grid[1], grid[2],
              pointData = {'velocity': reorderedU,
                           'pressure': rho})
