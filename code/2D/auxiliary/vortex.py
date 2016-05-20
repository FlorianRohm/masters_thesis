#!/usr/bin/python
from numpy import *

def rankinVortex(x,y,radius):
    r = sqrt(x**2+y**2)
    xVal = -y/radius
    yVal = x/radius
    if (r>radius):
        xVal = xVal*(radius/r)**2
        yVal = yVal*(radius/r)**2
    return [xVal,yVal]
