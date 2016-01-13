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
