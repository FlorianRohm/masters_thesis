from numpy import *

def centralMomentsFromDistributions_min (u, fin):
    ux = u[0, :, :]
    uy = u[1, :, :]
    (nx,ny) = ux.shape
    f = zeros((3,3,nx,ny))

    f[1,1] = fin[0]
    f[2,1] = fin[1]
    f[1,2] = fin[2]
    f[0,1] = fin[3]
    f[1,0] = fin[4]
    f[2,2] = fin[5]
    f[0,2] = fin[6]
    f[0,0] = fin[7]
    f[2,0] = fin[8]
    c_0_0 = f[0,0]              + f[0,1]           + f[0,2]
    c_0_1 = f[0,0]*(-1 - uy)    - f[0,1]*uy        + f[0,2]*(1 - uy)
    c_0_2 = f[0,0]*(-1 - uy)**2 + f[0,1]*uy**2     + f[0,2]*(1 - uy)**2
    c_1_0 = f[1,0]              + f[1,1]           + f[1,2]
    c_1_1 = f[1,0]*(-1 - uy)    - f[1,1]*uy        + f[1,2]*(1 - uy)
    c_1_2 = f[1,0]*(-1 - uy)**2 + f[1,1]*uy**2     + f[1,2]*(1 - uy)**2
    c_2_0 = f[2,0]              + f[2,1]           + f[2,2]
    c_2_1 = f[2,0]*(-1 - uy)    - f[2,1]*uy        + f[2,2]*(1 - uy)
    c_2_2 = f[2,0]*(-1 - uy)**2 + f[2,1]*uy**2     + f[2,2]*(1 - uy)**2

    c_20 = c_0_0*(-1 - ux)**2 + c_1_0*ux**2     + c_2_0*(1 - ux)**2
    c_11 = c_0_1*(-1 - ux)    - c_1_1*ux        + c_2_1*(1 - ux)
    c_02 = c_0_2              + c_1_2           + c_2_2

    return (c_11, c_20, c_02)
