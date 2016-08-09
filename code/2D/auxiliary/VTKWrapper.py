#!/usr/bin/python
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

from pyevtk.hl import gridToVTK
from numpy import *


def saveToVTK(velocity, rho, prefix, saveNumber, grid):
    name = "./" + prefix + "." + saveNumber
    gridToVTK(name, grid[0], grid[1], grid[2],
              pointData = {'velocity': velocity,
                           'pressure': rho})
