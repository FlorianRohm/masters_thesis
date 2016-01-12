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
            u[0] = c[2, i, j, 0] * f[2, i, j] + c[0, i, j, 0] * f[0, i, j]
            u[1] = c[i, 2, j, 1] * f[i, 2, j] + c[i, 0, j, 1] * f[i, 0, j]
            u[2] = c[i, j, 2, 2] * f[i, j, 2] + c[i, j, 0, 2] * f[i, j, 0]
    u = 3.0*u/rho
    return (rho, u)
