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

from numpy import *
from random import triangular, seed
seed(1503)

def rankinVortex(xDist,yDist,radius):
    r = sqrt(xDist**2+yDist**2)
    xVal = -yDist/radius
    yVal = xDist/radius
    if (r>radius):
        xVal = xVal*(radius/r)**2
        yVal = yVal*(radius/r)**2
    return [xVal,yVal]

def buildVortexField(nx,ny, nrOfPairs,halfDiameter):
    vel = zeros((2, nx, ny))
    for i in range(nrOfPairs):
        for j in range(nrOfPairs):
            c1x = (i*4 + 1)* halfDiameter*triangular(0.95,1.05)
            c2x = (i*4 + 3)* halfDiameter*triangular(0.95,1.05)
            c1y = (j*4 + 1)* halfDiameter*triangular(0.95,1.05)
            c2y = (j*4 + 3)* halfDiameter*triangular(0.95,1.05)
            #print("x1: {0}, y1: {1}, x2: {2}, y2: {3}".format(c1x, c2x, c1y, c2y))
            for x in range(nx):
                for y in range(ny):
                    for xoff in range(3):
                        xoff = xoff-1;
                        for yoff in range(3):
                            yoff = yoff-1;
                            vel[:,x,y] = vel[:,x,y] \
                            + rankinVortex(x-c1x + xoff*nx, y-c1y+ yoff*ny, halfDiameter/2.) \
                            + rankinVortex(x-c2x + xoff*nx, y-c2y+ yoff*ny, halfDiameter/2.) \
                            - rankinVortex(x-c2x + xoff*nx, y-c1y+ yoff*ny, halfDiameter/2.) \
                            - rankinVortex(x-c1x + xoff*nx, y-c2y+ yoff*ny, halfDiameter/2.)
    return vel
