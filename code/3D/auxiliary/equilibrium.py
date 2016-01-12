from numpy import *

q = 27

# lattice weights
t = ones((3,3,3))/54. # edge centers
t[1,1,1] = 8./27. # center
t[0,1,1] = t[1,0,1] = t[1,1,0] = \
t[2,1,1] = t[1,2,1] = t[1,1,2] = 2./27.  # face centers

t[0,0,0] = t[0,0,2] = t[0,2,0] = t[0,2,2] = \
t[2,0,0] = t[2,0,2] = t[2,2,0] = t[2,2,2] = 1./216. # corners

# Equilibrium distribution function.
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
