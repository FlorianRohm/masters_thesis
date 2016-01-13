from numpy import *
from equilibrium import equilibrium
from transformation.transformations import cumulantsFromDistribution, distributionFromCumulants
from getMacroValues import getMacroValues


def BGKCollide(fin, rho, u, omega):
    feq = equilibrium(rho, u)
    return fin - omega * (fin - feq)

def cumulantCollide(f, rho, u, omega):
    omega2  = 1
    omega3  = 1
    omega4  = 1
    omega5  = 1
    omega6  = 1
    omega7  = 1
    omega8  = 1
    omega9  = 1
    omega10 = 1

    rhoInv = 1./rho
    ux = u[0]
    uy = u[1]
    uz = u[2]

    ( C_000, C_001, C_002, C_010, C_011, C_012, C_020, C_021, C_022,
      C_100, C_101, C_102, C_110, C_111, C_112, C_120, C_121, C_122,
      C_200, C_201, C_202, C_210, C_211, C_212, C_220, C_221, C_222) = cumulantsFromDistribution(rhoInv, u, f)

    fTest = distributionFromCumulants(rhoInv, u, C_000, C_001, C_002, C_010, C_011, C_012, C_020, C_021, C_022,
                                                  C_100, C_101, C_102, C_110, C_111, C_112, C_120, C_121, C_122,
                                              C_200, C_201, C_202, C_210, C_211, C_212, C_220, C_221, C_222)

    (rho2, u2) = getMacroValues(fTest)
    ( C2_000, C2_001, C2_002, C2_010, C2_011, C2_012, C2_020, C2_021, C2_022,
      C2_100, C2_101, C2_102, C2_110, C2_111, C2_112, C2_120, C2_121, C2_122,
      C2_200, C2_201, C2_202, C2_210, C2_211, C2_212, C2_220, C2_221, C2_222) = cumulantsFromDistribution(1./rho2, u2, fTest)

    C_000_diff = C_000 - C2_000
    C_001_diff = C_001 - C2_001
    C_002_diff = C_002 - C2_002
    C_010_diff = C_010 - C2_010
    C_011_diff = C_011 - C2_011
    C_012_diff = C_012 - C2_012
    C_020_diff = C_020 - C2_020
    C_021_diff = C_021 - C2_021
    C_022_diff = C_022 - C2_022
    C_100_diff = C_100 - C2_100
    C_101_diff = C_101 - C2_101
    C_102_diff = C_102 - C2_102
    C_110_diff = C_110 - C2_110
    C_111_diff = C_111 - C2_111
    C_112_diff = C_112 - C2_112
    C_120_diff = C_120 - C2_120
    C_121_diff = C_121 - C2_121
    C_122_diff = C_122 - C2_122
    C_200_diff = C_200 - C2_200
    C_201_diff = C_201 - C2_201
    C_202_diff = C_202 - C2_202
    C_210_diff = C_210 - C2_210
    C_211_diff = C_211 - C2_211
    C_212_diff = C_212 - C2_212
    C_220_diff = C_220 - C2_220
    C_221_diff = C_221 - C2_221
    C_222_diff = C_222 - C2_222

    print "C_000 diff: {}".format(amax(C_000_diff))
    print "C_001 diff: {}".format(amax(C_001_diff))
    print "C_002 diff: {}".format(amax(C_002_diff))
    print "C_010 diff: {}".format(amax(C_010_diff))
    print "C_011 diff: {}".format(amax(C_011_diff))
    print "C_012 diff: {}".format(amax(C_012_diff))
    print "C_020 diff: {}".format(amax(C_020_diff))
    print "C_021 diff: {}".format(amax(C_021_diff))
    print "C_022 diff: {}".format(amax(C_022_diff))
    print "C_100 diff: {}".format(amax(C_100_diff))
    print "C_101 diff: {}".format(amax(C_101_diff))
    print "C_102 diff: {}".format(amax(C_102_diff))
    print "C_110 diff: {}".format(amax(C_110_diff))
    print "C_111 diff: {}".format(amax(C_111_diff))
    print "C_112 diff: {}".format(amax(C_112_diff))
    print "C_120 diff: {}".format(amax(C_120_diff))
    print "C_121 diff: {}".format(amax(C_121_diff))
    print "C_122 diff: {}".format(amax(C_122_diff))
    print "C_200 diff: {}".format(amax(C_200_diff))
    print "C_201 diff: {}".format(amax(C_201_diff))
    print "C_202 diff: {}".format(amax(C_202_diff))
    print "C_210 diff: {}".format(amax(C_210_diff))
    print "C_211 diff: {}".format(amax(C_211_diff))
    print "C_212 diff: {}".format(amax(C_212_diff))
    print "C_220 diff: {}".format(amax(C_220_diff))
    print "C_221 diff: {}".format(amax(C_221_diff))
    print "C_222 diff: {}".format(amax(C_222_diff))

    rhodiff = rho-rho2
    udiff = u-u2

    print "rho: {}\nu:   {}\n".format(amax(rhodiff), amax(udiff))

    CS_000 = C_000
    CS_100 = C_100
    CS_010 = C_010
    CS_001 = C_001

    CS_110 = (1-omega)*C_110
    CS_101 = (1-omega)*C_101
    CS_011 = (1-omega)*C_011

    # (C_000 =) k_000 in the paper, not sure, why he wrote k
    Dxux = - 0.5*omega*rhoInv*(2*C_200 - C_020 - C_002) - 0.5*1*rhoInv*(C_200 + C_020 + C_002 - C_000)
    Dyuy = Dxux + 1.5*omega*rhoInv*(C_200 - C_020)
    Dzuz = Dxux + 1.5*omega*rhoInv*(C_200 - C_002)

    CS_200__m__CS020 = (1 - omega)*(C_200 - C_020) - 3*rho*(1 - 0.5*omega)*(Dxux*ux**2 - Dyuy*uy**2)
    CS_200__m__CS002 = (1 - omega)*(C_200 - C_002) - 3*rho*(1 - 0.5*omega)*(Dxux*ux**2 - Dzuz*uz**2)
    CS_200__p__CS020__p__CS_002 = omega2*C_000 + (1-omega2)*(C_200 + C_020 + C_002) \
                                  - 3*rho*(1 - 0.5*omega2)*(Dxux*ux**2 + Dyuy*uy**2 + Dzuz*uz**2)

    CS_200 = (CS_200__m__CS020 + CS_200__m__CS002 + CS_200__p__CS020__p__CS_002)/3.
    CS_020 = CS_200 - CS_200__m__CS020
    CS_002 = CS_200 - CS_200__m__CS002


    CS_120__p__CS_102 = (1 - omega3)*(C_120 + C_102)
    CS_210__p__CS_012 = (1 - omega3)*(C_210 + C_012)
    CS_201__p__CS_021 = (1 - omega3)*(C_201 + C_021)
    CS_120__m__CS_102 = (1 - omega4)*(C_120 - C_102)
    CS_210__m__CS_012 = (1 - omega4)*(C_210 - C_012)
    CS_201__m__CS_021 = (1 - omega4)*(C_201 - C_021)

    CS_120 = 0.5*(CS_120__p__CS_102 + CS_120__m__CS_102)
    CS_102 = 0.5*(CS_120__p__CS_102 - CS_120__m__CS_102)
    CS_012 = 0.5*(CS_210__p__CS_012 - CS_210__m__CS_012)
    CS_210 = 0.5*(CS_210__p__CS_012 + CS_210__m__CS_012)
    CS_201 = 0.5*(CS_201__p__CS_021 + CS_201__m__CS_021)
    CS_021 = 0.5*(CS_201__p__CS_021 - CS_201__m__CS_021)

    CS_111 = (1 - omega5) * C_111


    CS_220__m__2CS_202__p__CS_022 = (1 - omega6)*(C_220 - 2*C_202 + C_022)
    CS_220__p__CS_202__m__2CS_022 = (1 - omega6)*(C_220 + C_202 - 2*C_022)
    CS_220__p__CS_202__p__CS_022  = (1 - omega7)*(C_220 + C_202 + C_022)

    CS_220 = (CS_220__m__2CS_202__p__CS_022 + CS_220__p__CS_202__m__2CS_022 + CS_220__p__CS_202__p__CS_022)/3.
    CS_202 = (CS_220__p__CS_202__p__CS_022 - CS_220__m__2CS_202__p__CS_022)/3.
    CS_022 = (CS_220__p__CS_202__p__CS_022 - CS_220__p__CS_202__m__2CS_022)/3.

    CS_211 = (1 - omega8 )* C_211
    CS_121 = (1 - omega8 )* C_121
    CS_112 = (1 - omega8 )* C_112
    CS_221 = (1 - omega9 )* C_221
    CS_212 = (1 - omega9 )* C_212
    CS_122 = (1 - omega9 )* C_122
    CS_222 = (1 - omega10 )* C_222

    return distributionFromCumulants(rhoInv, u, CS_000, CS_001, CS_002, CS_010, CS_011, CS_012, CS_020, CS_021, CS_022,
                                                CS_100, CS_101, CS_102, CS_110, CS_111, CS_112, CS_120, CS_121, CS_122,
                                                CS_200, CS_201, CS_202, CS_210, CS_211, CS_212, CS_220, CS_221, CS_222)
