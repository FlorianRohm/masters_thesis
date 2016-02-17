from numpy import *
from transformations.transforms import *

def BGKCollide(fin, feq, omega):
    return fin - omega * (fin - feq)

def cumulantCollide(fin, rho, u, omega):
    (K_00, K_10, K_01, K_11, K_20, K_02, K_21, K_12, K_22) = normalizedCumulantsFromDistributions (rho, u, fin)

    #print "\ncollide"
    # collision
    K_00_p = K_00
    K_10_p = K_10
    K_01_p = K_01

    K_11_p = (1-omega)*K_11

    K_20_p = 1/3 + 0.5*(1-omega)*(K_20 - K_02)
    K_02_p = 1/3 - 0.5*(1-omega)*(K_20 - K_02)

    K_21_p = 0
    K_12_p = 0

    K_22_p = 0

    return distributionsFromNormalizedCumulants (rho, u, K_00_p, K_10_p, K_01_p, K_11_p, K_20_p, K_02_p, K_21_p, K_12_p, K_22_p)

def cumulantCollideAll(fin, rho, u, omega1, omega2, omega3, omega4):
    (K_00, K_10, K_01, K_11, K_20, K_02, K_21, K_12, K_22) = normalizedCumulantsFromDistributions (rho, u, fin)

    # collision
    K_00_p = K_00
    K_10_p = K_10
    K_01_p = K_01

    K_11_p = (1-omega1)*K_11

    K_20_p = 0.5 * (K_20 + K_02 + omega2*(2/3 - K_20 - K_02) + (1-omega1)*(K_20 - K_02))
    K_02_p = 0.5 * (K_20 + K_02 + omega2*(2/3 - K_20 - K_02) - (1-omega1)*(K_20 - K_02))

    K_21_p = (1-omega3)*K_21
    K_12_p = (1-omega3)*K_12

    K_22_p = (1-omega4)*K_22

    return distributionsFromNormalizedCumulants (rho, u, K_00_p, K_10_p, K_01_p, K_11_p, K_20_p, K_02_p, K_21_p, K_12_p, K_22_p)


def centralMomentSRT(fin, feq, u, omega):
    # central moments
    (c_00, c_10, c_01, c_11, c_20, c_02, c_21, c_12, c_22) = centralMomentsFromDistributions (u, fin)

    # transform the equilibrium function
    (c_00_eq, c_10_eq, c_01_eq, c_11_eq, c_20_eq, c_02_eq, c_21_eq, c_12_eq, c_22_eq) = centralMomentsFromDistributions (u, feq)

    # collision
    c_00_p = c_00 + omega*(c_00_eq - c_00)
    c_10_p = c_10 + omega*(c_10_eq - c_10)
    c_01_p = c_01 + omega*(c_01_eq - c_01)

    c_11_p = c_11 + omega*(c_11_eq - c_11)

    c_20_p = c_20 + omega*(c_20_eq - c_20)
    c_02_p = c_02 + omega*(c_02_eq - c_02)

    c_21_p = c_21 + omega*(c_21_eq - c_21)
    c_12_p = c_12 + omega*(c_12_eq - c_12)

    c_22_p = c_22 + omega*(c_22_eq - c_22)

    return distributionsFromCentralMoments (u, c_00_p, c_10_p, c_01_p, c_11_p, c_20_p, c_02_p, c_21_p, c_12_p, c_22_p)
