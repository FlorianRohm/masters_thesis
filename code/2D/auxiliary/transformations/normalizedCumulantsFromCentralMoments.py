from numpy import *

def normalizedCumulantsFromCentralMoments (rho, u, c_00, c_10, c_01, c_11, c_20, c_02, c_21, c_12, c_22):

    ux = u[0, :, :]
    uy = u[1, :, :]

    K_22 = c_22 - 2*c_11*c_11/rho - c_20*c_02/rho

    return (log(rho)*rho, ux*rho, uy*rho, c_11, c_20, c_02, c_21, c_12, K_22)
