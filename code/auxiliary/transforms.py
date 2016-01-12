from numpy import *

def centralMomentFromDistribution_Long(u, fin):
    ux = u[0, :, :]
    uy = u[1, :, :]

    fin_00 = fin[0, :, :]
    fin_10 = fin[1, :, :]
    fin_01 = fin[2, :, :]
    fin_n10 = fin[3, :, :]
    fin_0n1 = fin[4, :, :]
    fin_11 = fin[5, :, :]
    fin_n11 = fin[6, :, :]
    fin_n1n1 = fin[7, :, :]
    fin_1n1 = fin[8, :, :]

    c_00 = \
    fin_n1n1 + \
    fin_n10 + \
    fin_n11 + \
    fin_0n1 + \
    fin_00 + \
    fin_01 + \
    fin_1n1 + \
    fin_10 + \
    fin_11

    c_10 = \
    (-1 - ux) * fin_n1n1 + \
    (-1 - ux) * fin_n10 + \
    (-1 - ux) * fin_n11 + \
    (- ux)    * fin_0n1 + \
    (- ux)    * fin_00 + \
    (- ux)    * fin_01 + \
    (1 - ux)  * fin_1n1 + \
    (1 - ux)  * fin_10 + \
    (1 - ux)  * fin_11

    c_01 = \
    (-1 - uy) * fin_n1n1 + \
    ( - uy)   * fin_n10 + \
    (1 - uy)  * fin_n11 + \
    (-1 - uy) * fin_0n1 + \
    ( - uy)   * fin_00 + \
    (1 - uy)  * fin_01 + \
    (-1 - uy) * fin_1n1 + \
    ( - uy)   * fin_10 + \
    (1 - uy)  * fin_11

    c_11 = \
    (-1 - ux) * (-1 - uy) * fin_n1n1 + \
    (-1 - ux) * ( - uy)   * fin_n10 + \
    (-1 - ux) * (1 - uy)  * fin_n11 + \
    (- ux)    * (-1 - uy) * fin_0n1 + \
    (- ux)    * ( - uy)   * fin_00 + \
    (- ux)    * (1 - uy)  * fin_01 + \
    (1 - ux)  * (-1 - uy) * fin_1n1 + \
    (1 - ux)  * ( - uy)   * fin_10 + \
    (1 - ux)  * (1 - uy)  * fin_11

    c_20 = \
    square((-1 - ux)) * fin_n1n1 + \
    square((-1 - ux)) * fin_n10 + \
    square((-1 - ux)) * fin_n11 + \
    square((- ux))    * fin_0n1 + \
    square((- ux))    * fin_00 + \
    square((- ux))    * fin_01 + \
    square((1 - ux))  * fin_1n1 + \
    square((1 - ux))  * fin_10 + \
    square((1 - ux))  * fin_11

    c_02 = \
    square((-1 - uy)) * fin_n1n1 + \
    square(( - uy))   * fin_n10 + \
    square((1 - uy))  * fin_n11 + \
    square((-1 - uy)) * fin_0n1 + \
    square(( - uy))   * fin_00 + \
    square((1 - uy))  * fin_01 + \
    square((-1 - uy)) * fin_1n1 + \
    square(( - uy))   * fin_10 + \
    square((1 - uy))  * fin_11

    c_21 = \
    square((-1 - ux)) * (-1 - uy) * fin_n1n1 + \
    square((-1 - ux)) * ( - uy)   * fin_n10 + \
    square((-1 - ux)) * (1 - uy)  * fin_n11 + \
    square((- ux))    * (-1 - uy) * fin_0n1 + \
    square((- ux))    * ( - uy)   * fin_00 + \
    square((- ux))    * (1 - uy)  * fin_01 + \
    square((1 - ux))  * (-1 - uy) * fin_1n1 + \
    square((1 - ux))  * ( - uy)   * fin_10 + \
    square((1 - ux))  * (1 - uy)  * fin_11

    c_12 = \
    (-1 - ux) * square((-1 - uy)) * fin_n1n1 + \
    (-1 - ux) * square(( - uy))   * fin_n10 + \
    (-1 - ux) * square((1 - uy))  * fin_n11 + \
    (- ux)    * square((-1 - uy)) * fin_0n1 + \
    (- ux)    * square(( - uy))   * fin_00 + \
    (- ux)    * square((1 - uy))  * fin_01 + \
    (1 - ux)  * square((-1 - uy)) * fin_1n1 + \
    (1 - ux)  * square(( - uy))   * fin_10 + \
    (1 - ux)  * square((1 - uy))  * fin_11

    c_22 = \
    square((-1 - ux)) * square((-1 - uy)) * fin_n1n1 + \
    square((-1 - ux)) * square(( - uy))   * fin_n10 + \
    square((-1 - ux)) * square((1 - uy))  * fin_n11 + \
    square((- ux))    * square((-1 - uy)) * fin_0n1 + \
    square((- ux))    * square(( - uy))   * fin_00 + \
    square((- ux))    * square((1 - uy))  * fin_01 + \
    square((1 - ux))  * square((-1 - uy)) * fin_1n1 + \
    square((1 - ux))  * square(( - uy))   * fin_10 + \
    square((1 - ux))  * square((1 - uy))  * fin_11

    return (c_00, c_10, c_01, c_11, c_20, c_02, c_21, c_12, c_22)


def centralMomentFromDistribution (u, fin):
    ux = u[0, :, :]
    uy = u[1, :, :]

    # c_{i beta}
    # c_{-1 beta}
    cDash_n1_0 =       fin[3, :, :] +                 fin[7, :, :] +               fin[6, :, :]
    cDash_n1_1 =   -uy*fin[3, :, :] +         (-1-uy)*fin[7, :, :] +        (1-uy)*fin[6, :, :]
    cDash_n1_2 = uy*uy*fin[3, :, :] + (-1-uy)*(-1-uy)*fin[7, :, :] + (1-uy)*(1-uy)*fin[6, :, :]

    # c_{0 beta}
    cDash_0_0 =       fin[0, :, :] +                 fin[4, :, :] +               fin[2, :, :]
    cDash_0_1 =   -uy*fin[0, :, :] +         (-1-uy)*fin[4, :, :] +        (1-uy)*fin[2, :, :]
    cDash_0_2 = uy*uy*fin[0, :, :] + (-1-uy)*(-1-uy)*fin[4, :, :] + (1-uy)*(1-uy)*fin[2, :, :]

    # c_{1 beta}
    cDash_1_0 =       fin[1, :, :] +                 fin[8, :, :] +               fin[5, :, :]
    cDash_1_1 =   -uy*fin[1, :, :] +         (-1-uy)*fin[8, :, :] +        (1-uy)*fin[5, :, :]
    cDash_1_2 = uy*uy*fin[1, :, :] + (-1-uy)*(-1-uy)*fin[8, :, :] + (1-uy)*(1-uy)*fin[5, :, :]

    # c{alpha beta}
    # c{alpha 0}
    c_00 =                 cDash_n1_0 +       cDash_0_0 +               cDash_1_0
    c_10 =         (-1-ux)*cDash_n1_0 -    ux*cDash_0_0 +        (1-ux)*cDash_1_0
    c_20 = (-1-ux)*(-1-ux)*cDash_n1_0 + ux*ux*cDash_0_0 + (1-ux)*(1-ux)*cDash_1_0

    # c{alpha 1}
    c_01 =                 cDash_n1_1 +       cDash_0_1 +               cDash_1_1
    c_11 =         (-1-ux)*cDash_n1_1 -    ux*cDash_0_1 +        (1-ux)*cDash_1_1
    c_21 = (-1-ux)*(-1-ux)*cDash_n1_1 + ux*ux*cDash_0_1 + (1-ux)*(1-ux)*cDash_1_1

    # c{alpha 2}
    c_02 =                 cDash_n1_2 +       cDash_0_2 +               cDash_1_2
    c_12 =         (-1-ux)*cDash_n1_2 -    ux*cDash_0_2 +        (1-ux)*cDash_1_2
    c_22 = (-1-ux)*(-1-ux)*cDash_n1_2 + ux*ux*cDash_0_2 + (1-ux)*(1-ux)*cDash_1_2
    return (c_00, c_10, c_01, c_11, c_20, c_02, c_21, c_12, c_22)

def distributionsFromCentralMoments (u, c_00, c_10, c_01, c_11, c_20, c_02, c_21, c_12, c_22):
    ux = u[0, :, :]
    uy = u[1, :, :]

    # backward transformation
    # cDash_n1_0 = -0.5*(ux*(1 - ux))*c_00 - 0.5*(1-2*ux)*c_10 + 0.5*c_20
    # cDash_0_0  =          (1-ux*ux)*c_00 -         2*ux*c_10 -     c_20
    # cDash_1_0  =  0.5*(ux*(1 + ux))*c_00 + 0.5*(1+2*ux)*c_10 + 0.5*c_20

    cDash_n1_0 = -0.5*(ux*(1 - ux))*c_00 - 0.5*(1-2*ux)*c_10 + 0.5*c_20
    cDash_0_0  =          (1-ux*ux)*c_00 -         2*ux*c_10 -     c_20
    cDash_1_0  =  0.5*(ux*(1 + ux))*c_00 + 0.5*(1+2*ux)*c_10 + 0.5*c_20

    cDash_n1_1 = -0.5*(ux*(1 - ux))*c_01 - 0.5*(1-2*ux)*c_11 + 0.5*c_21
    cDash_0_1  =          (1-ux*ux)*c_01 -         2*ux*c_11 -     c_21
    cDash_1_1  =  0.5*(ux*(1 + ux))*c_01 + 0.5*(1+2*ux)*c_11 + 0.5*c_21

    cDash_n1_2 = -0.5*(ux*(1 - ux))*c_02 - 0.5*(1-2*ux)*c_12 + 0.5*c_22
    cDash_0_2  =          (1-ux*ux)*c_02 -         2*ux*c_12 -     c_22
    cDash_1_2  =  0.5*(ux*(1 + ux))*c_02 + 0.5*(1+2*ux)*c_12 + 0.5*c_22

    # post collision distributions
    # f__i_n1 = -0.5*(uy*(1 - uy))*cDash_i_0 - 0.5*(1-2*uy)*cDash_i_1 + 0.5*cDash_i_2
    # f__i_0  =        (1 - uy*uy)*cDash_i_0 -         2*uy*cDash_i_1 -     cDash_i_2
    # f__i_1  = -0.5*(uy*(1 + uy))*cDash_i_0 + 0.5*(1+2*uy)*cDash_i_1 + 0.5*cDash_i_2

    f__n1_n1 = -0.5*(uy*(1 - uy))*cDash_n1_0 - 0.5*(1-2*uy)*cDash_n1_1 + 0.5*cDash_n1_2
    f__n1_0  =        (1 - uy*uy)*cDash_n1_0 -         2*uy*cDash_n1_1 -     cDash_n1_2
    f__n1_1  = -0.5*(uy*(1 + uy))*cDash_n1_0 + 0.5*(1+2*uy)*cDash_n1_1 + 0.5*cDash_n1_2

    f__0_n1 = -0.5*(uy*(1 - uy))*cDash_0_0 - 0.5*(1-2*uy)*cDash_0_1 + 0.5*cDash_0_2
    f__0_0  =        (1 - uy*uy)*cDash_0_0 -         2*uy*cDash_0_1 -     cDash_0_2
    f__0_1  = -0.5*(uy*(1 + uy))*cDash_0_0 + 0.5*(1+2*uy)*cDash_0_1 + 0.5*cDash_0_2

    f__1_n1 = -0.5*(uy*(1 - uy))*cDash_1_0 - 0.5*(1-2*uy)*cDash_1_1 + 0.5*cDash_1_2
    f__1_0  =        (1 - uy*uy)*cDash_1_0 -         2*uy*cDash_1_1 -     cDash_1_2
    f__1_1  = -0.5*(uy*(1 + uy))*cDash_1_0 + 0.5*(1+2*uy)*cDash_1_1 + 0.5*cDash_1_2

    return array([f__0_0, f__1_0, f__0_1, f__n1_0, f__0_n1, f__1_1, f__n1_1, f__n1_n1, f__1_n1])

def normalizedCumulantsFromCentralMoments (rho, c_00, c_10, c_01, c_11, c_20, c_02, c_21, c_12, c_22):
    # only K_22 differs from the central moments
    K_22 = c_22 - 2*c_11*c_11/rho - c_20*c_02/rho

    return (c_00, c_10, c_01, c_11, c_20, c_02, c_21, c_12, K_22)

def centralMomentsFromNormalizedCumulants (rho, K_00, K_10, K_01, K_11, K_20, K_02, K_21, K_12, K_22):
    # Transformation to central moments
    c_22 = K_22 + 2*K_11*K_11/rho + K_20*K_02/rho

    return (K_00, K_10, K_01, K_11, K_20, K_02, K_21, K_12, c_22)

def normalizedCumulantsFromDistributions(rho, u, fin):
    (c_00, c_10, c_01, c_11, c_20, c_02, c_21, c_12, c_22) = centralMomentFromDistribution (u, fin)

    return normalizedCumulantsFromCentralMoments (rho,c_00, c_10, c_01, c_11, c_20, c_02, c_21, c_12, c_22)

def normalizedCumulantsFromDistributions_Long(rho, u, fin):
    (c_00, c_10, c_01, c_11, c_20, c_02, c_21, c_12, c_22) = centralMomentFromDistribution_Long (u, fin)

    return normalizedCumulantsFromCentralMoments (rho,c_00, c_10, c_01, c_11, c_20, c_02, c_21, c_12, c_22)

def distributionsFromNormalizedCumulants(rho, u, K_00, K_10, K_01, K_11, K_20, K_02, K_21, K_12, K_22):
    (c_00, c_10, c_01, c_11, c_20, c_02, c_21, c_12, c_22) = \
    normalizedCumulantsFromCentralMoments (rho, K_00, K_10, K_01, K_11, K_20, K_02, K_21, K_12, K_22)

    return distributionsFromCentralMoments (u, c_00, c_10, c_01, c_11, c_20, c_02, c_21, c_12, c_22)
