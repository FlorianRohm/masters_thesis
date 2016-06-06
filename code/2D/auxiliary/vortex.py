#!/usr/bin/python
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
