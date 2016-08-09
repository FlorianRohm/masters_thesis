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

from scipy import misc
import numpy


def importImage(path):
    image = misc.imread(path)
    image = image[:, :, 0]
    image = image < 50.
    image = numpy.swapaxes(image, 1, 0)
    return image


def flipY(image):
    return numpy.fliplr(image)


def setSubmatrixAt(matrix, submatrix, startX, startY):
    dx = submatrix.shape[0]
    dy = submatrix.shape[1]
    matrix[startX: startX + dx, startY: startY+dy] = submatrix
    return matrix

if __name__ == "__main__":
    import sys
    importImage(sys.argv[1])
