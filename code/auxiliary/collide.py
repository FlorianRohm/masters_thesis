from numpy import *

def BGKCollide(fin, feq, omega):
    return fin - omega * (fin - feq)


def cumulantCollide(fin, u, rho, omega):
    # central moments
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

    # only K22 differs from the central moments
    K_22 = c_22 - 2*c_11*c_11/rho - c_20*c_02/rho

    #print "\ncollide"
    # collision
    c_00_post = c_00
    c_10_post = c_10
    c_01_post = c_01

    c_11_post = (1-omega)*c_11

    c_20_post = 1/3 + 0.5*(1-omega)*(c_20 - c_02)
    c_02_post = 1/3 - 0.5*(1-omega)*(c_20 - c_02)

    c_21_post = 0
    c_12_post = 0

    K_22_post = 0

    # no Collision
    # c_00_post = c_00
    # c_10_post = c_10
    # c_01_post = c_01
    # c_11_post = c_11
    # c_20_post = c_20
    # c_02_post = c_02

    # c_21_post = c_21
    # c_12_post = c_12

    # K_22_post = K_22

    # Transformation to central moments
    c_22_post = K_22_post + 2*c_11_post*c_11_post/rho + c_20_post*c_02_post/rho


    # backward transformation
    # cDash_n1_0_post = -0.5*(ux*(1 - ux))*c_00_post - 0.5*(1-2*ux)*c_10_post + 0.5*c_20_post
    # cDash_0_0_post  =          (1-ux*ux)*c_00_post -         2*ux*c_10_post -     c_20_post
    # cDash_1_0_post  =  0.5*(ux*(1 + ux))*c_00_post + 0.5*(1+2*ux)*c_10_post + 0.5*c_20_post

    cDash_n1_0_post = -0.5*(ux*(1 - ux))*c_00_post - 0.5*(1-2*ux)*c_10_post + 0.5*c_20_post
    cDash_0_0_post  =          (1-ux*ux)*c_00_post -         2*ux*c_10_post -     c_20_post
    cDash_1_0_post  =  0.5*(ux*(1 + ux))*c_00_post + 0.5*(1+2*ux)*c_10_post + 0.5*c_20_post

    cDash_n1_1_post = -0.5*(ux*(1 - ux))*c_01_post - 0.5*(1-2*ux)*c_11_post + 0.5*c_21_post
    cDash_0_1_post  =          (1-ux*ux)*c_01_post -         2*ux*c_11_post -     c_21_post
    cDash_1_1_post  =  0.5*(ux*(1 + ux))*c_01_post + 0.5*(1+2*ux)*c_11_post + 0.5*c_21_post

    cDash_n1_2_post = -0.5*(ux*(1 - ux))*c_02_post - 0.5*(1-2*ux)*c_12_post + 0.5*c_22_post
    cDash_0_2_post  =          (1-ux*ux)*c_02_post -         2*ux*c_12_post -     c_22_post
    cDash_1_2_post  =  0.5*(ux*(1 + ux))*c_02_post + 0.5*(1+2*ux)*c_12_post + 0.5*c_22_post

    # post collision distributions
    # f__i_n1_post = -0.5*(uy*(1 - uy))*cDash_i_0_post - 0.5*(1-2*uy)*cDash_i_1_post + 0.5*cDash_i_2_post
    # f__i_0_post  =        (1 - uy*uy)*cDash_i_0_post -         2*uy*cDash_i_1_post -     cDash_i_2_post
    # f__i_1_post  = -0.5*(uy*(1 + uy))*cDash_i_0_post + 0.5*(1+2*uy)*cDash_i_1_post + 0.5*cDash_i_2_post

    f__n1_n1_post = -0.5*(uy*(1 - uy))*cDash_n1_0_post - 0.5*(1-2*uy)*cDash_n1_1_post + 0.5*cDash_n1_2_post
    f__n1_0_post  =        (1 - uy*uy)*cDash_n1_0_post -         2*uy*cDash_n1_1_post -     cDash_n1_2_post
    f__n1_1_post  = -0.5*(uy*(1 + uy))*cDash_n1_0_post + 0.5*(1+2*uy)*cDash_n1_1_post + 0.5*cDash_n1_2_post

    f__0_n1_post = -0.5*(uy*(1 - uy))*cDash_0_0_post - 0.5*(1-2*uy)*cDash_0_1_post + 0.5*cDash_0_2_post
    f__0_0_post  =        (1 - uy*uy)*cDash_0_0_post -         2*uy*cDash_0_1_post -     cDash_0_2_post
    f__0_1_post  = -0.5*(uy*(1 + uy))*cDash_0_0_post + 0.5*(1+2*uy)*cDash_0_1_post + 0.5*cDash_0_2_post

    f__1_n1_post = -0.5*(uy*(1 - uy))*cDash_1_0_post - 0.5*(1-2*uy)*cDash_1_1_post + 0.5*cDash_1_2_post
    f__1_0_post  =        (1 - uy*uy)*cDash_1_0_post -         2*uy*cDash_1_1_post -     cDash_1_2_post
    f__1_1_post  = -0.5*(uy*(1 + uy))*cDash_1_0_post + 0.5*(1+2*uy)*cDash_1_1_post + 0.5*cDash_1_2_post

    return array([f__0_0_post, f__1_0_post, f__0_1_post, f__n1_0_post, f__0_n1_post,   f__1_1_post, f__n1_1_post, f__n1_n1_post, f__1_n1_post])


def cumulantCollideAll(fin, u, rho, omega1, omega2, omega3, omega4):
    # central moments
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

    # only K22 differs from the central moments
    K_22 = c_22 - 2*c_11*c_11/rho - c_20*c_02/rho

    # collision
    c_00_post = c_00
    c_10_post = c_10
    c_01_post = c_01

    c_11_post = (1-omega1)*c_11

    c_20_post = 0.5 * (c_20 + c_02 + omega2*(2/3 - c_20 - c_02) + (1-omega1)*(c_20 - c_02))
    c_02_post = 0.5 * (c_20 + c_02 + omega2*(2/3 - c_20 - c_02) - (1-omega1)*(c_20 - c_02))

    c_21_post = (1-omega3)*c_21
    c_12_post = (1-omega3)*c_12

    K_22_post = (1-omega4)*K_22

    # Transformation to central moments
    c_22_post = K_22_post + 2*c_11_post*c_11_post/rho + c_20_post*c_02_post/rho

    # backward transformation
    cDash_n1_0_post = -0.5*(ux*(1 - ux))*c_00_post - 0.5*(1-2*ux)*c_10_post + 0.5*c_20_post
    cDash_0_0_post  =          (1-ux*ux)*c_00_post -         2*ux*c_10_post -     c_20_post
    cDash_1_0_post  =  0.5*(ux*(1 + ux))*c_00_post + 0.5*(1+2*ux)*c_10_post + 0.5*c_20_post

    cDash_n1_1_post = -0.5*(ux*(1 - ux))*c_01_post - 0.5*(1-2*ux)*c_11_post + 0.5*c_21_post
    cDash_0_1_post  =          (1-ux*ux)*c_01_post -         2*ux*c_11_post -     c_21_post
    cDash_1_1_post  =  0.5*(ux*(1 + ux))*c_01_post + 0.5*(1+2*ux)*c_11_post + 0.5*c_21_post

    cDash_n1_2_post = -0.5*(ux*(1 - ux))*c_02_post - 0.5*(1-2*ux)*c_12_post + 0.5*c_22_post
    cDash_0_2_post  =          (1-ux*ux)*c_02_post -         2*ux*c_12_post -     c_22_post
    cDash_1_2_post  =  0.5*(ux*(1 + ux))*c_02_post + 0.5*(1+2*ux)*c_12_post + 0.5*c_22_post

    # post collision distributions
    f__n1_n1_post = -0.5*(uy*(1 - uy))*cDash_n1_0_post - 0.5*(1-2*uy)*cDash_n1_1_post + 0.5*cDash_n1_2_post
    f__n1_0_post  =        (1 - uy*uy)*cDash_n1_0_post -         2*uy*cDash_n1_1_post -     cDash_n1_2_post
    f__n1_1_post  = -0.5*(uy*(1 + uy))*cDash_n1_0_post + 0.5*(1+2*uy)*cDash_n1_1_post + 0.5*cDash_n1_2_post

    f__0_n1_post = -0.5*(uy*(1 - uy))*cDash_0_0_post - 0.5*(1-2*uy)*cDash_0_1_post + 0.5*cDash_0_2_post
    f__0_0_post  =        (1 - uy*uy)*cDash_0_0_post -         2*uy*cDash_0_1_post -     cDash_0_2_post
    f__0_1_post  = -0.5*(uy*(1 + uy))*cDash_0_0_post + 0.5*(1+2*uy)*cDash_0_1_post + 0.5*cDash_0_2_post

    f__1_n1_post = -0.5*(uy*(1 - uy))*cDash_1_0_post - 0.5*(1-2*uy)*cDash_1_1_post + 0.5*cDash_1_2_post
    f__1_0_post  =        (1 - uy*uy)*cDash_1_0_post -         2*uy*cDash_1_1_post -     cDash_1_2_post
    f__1_1_post  = -0.5*(uy*(1 + uy))*cDash_1_0_post + 0.5*(1+2*uy)*cDash_1_1_post + 0.5*cDash_1_2_post

    return array([f__0_0_post, f__1_0_post, f__0_1_post, f__n1_0_post, f__0_n1_post,   f__1_1_post, f__n1_1_post, f__n1_n1_post, f__1_n1_post])



def centralMomentSRT(fin, feq, u, omega):
    # central moments
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

    # transform the equilibrium function
    # c_{i beta}
    # c_{-1 beta}
    cDash_n1_0_eq =       feq[3, :, :] +                 feq[7, :, :] +               feq[6, :, :]
    cDash_n1_1_eq =   -uy*feq[3, :, :] +         (-1-uy)*feq[7, :, :] +        (1-uy)*feq[6, :, :]
    cDash_n1_2_eq = uy*uy*feq[3, :, :] + (-1-uy)*(-1-uy)*feq[7, :, :] + (1-uy)*(1-uy)*feq[6, :, :]

    # c_{0 beta}
    cDash_0_0_eq =       feq[0, :, :] +                 feq[4, :, :] +               feq[2, :, :]
    cDash_0_1_eq =   -uy*feq[0, :, :] +         (-1-uy)*feq[4, :, :] +        (1-uy)*feq[2, :, :]
    cDash_0_2_eq = uy*uy*feq[0, :, :] + (-1-uy)*(-1-uy)*feq[4, :, :] + (1-uy)*(1-uy)*feq[2, :, :]

    # c_{1 beta}
    cDash_1_0_eq =       feq[1, :, :] +                 feq[8, :, :] +               feq[5, :, :]
    cDash_1_1_eq =   -uy*feq[1, :, :] +         (-1-uy)*feq[8, :, :] +        (1-uy)*feq[5, :, :]
    cDash_1_2_eq = uy*uy*feq[1, :, :] + (-1-uy)*(-1-uy)*feq[8, :, :] + (1-uy)*(1-uy)*feq[5, :, :]

    # c{alpha beta}
    # c{alpha 0}
    c_00_eq =                 cDash_n1_0_eq +       cDash_0_0_eq +               cDash_1_0_eq
    c_10_eq =         (-1-ux)*cDash_n1_0_eq -    ux*cDash_0_0_eq +        (1-ux)*cDash_1_0_eq
    c_20_eq = (-1-ux)*(-1-ux)*cDash_n1_0_eq + ux*ux*cDash_0_0_eq + (1-ux)*(1-ux)*cDash_1_0_eq

    # c{alpha 1}
    c_01_eq =                 cDash_n1_1_eq +       cDash_0_1_eq +               cDash_1_1_eq
    c_11_eq =         (-1-ux)*cDash_n1_1_eq -    ux*cDash_0_1_eq +        (1-ux)*cDash_1_1_eq
    c_21_eq = (-1-ux)*(-1-ux)*cDash_n1_1_eq + ux*ux*cDash_0_1_eq + (1-ux)*(1-ux)*cDash_1_1_eq

    # c{alpha 2}
    c_02_eq =                 cDash_n1_2_eq +       cDash_0_2_eq +               cDash_1_2_eq
    c_12_eq =         (-1-ux)*cDash_n1_2_eq -    ux*cDash_0_2_eq +        (1-ux)*cDash_1_2_eq
    c_22_eq = (-1-ux)*(-1-ux)*cDash_n1_2_eq + ux*ux*cDash_0_2_eq + (1-ux)*(1-ux)*cDash_1_2_eq


    # collision
    c_00_post = c_00 + omega*(c_00_eq - c_00)
    c_10_post = c_10 + omega*(c_10_eq - c_10)
    c_01_post = c_01 + omega*(c_01_eq - c_01)

    c_11_post = c_11 + omega*(c_11_eq - c_11)

    c_20_post = c_20 + omega*(c_20_eq - c_20)
    c_02_post = c_02 + omega*(c_02_eq - c_02)

    c_21_post = c_21 + omega*(c_21_eq - c_21)
    c_12_post = c_12 + omega*(c_12_eq - c_12)

    c_22_post = c_22 + omega*(c_22_eq - c_22)

    # backward transformation
    cDash_n1_0_post = -0.5*(ux*(1 - ux))*c_00_post - 0.5*(1-2*ux)*c_10_post + 0.5*c_20_post
    cDash_0_0_post  =          (1-ux*ux)*c_00_post -         2*ux*c_10_post -     c_20_post
    cDash_1_0_post  =  0.5*(ux*(1 + ux))*c_00_post + 0.5*(1+2*ux)*c_10_post + 0.5*c_20_post

    cDash_n1_1_post = -0.5*(ux*(1 - ux))*c_01_post - 0.5*(1-2*ux)*c_11_post + 0.5*c_21_post
    cDash_0_1_post  =          (1-ux*ux)*c_01_post -         2*ux*c_11_post -     c_21_post
    cDash_1_1_post  =  0.5*(ux*(1 + ux))*c_01_post + 0.5*(1+2*ux)*c_11_post + 0.5*c_21_post

    cDash_n1_2_post = -0.5*(ux*(1 - ux))*c_02_post - 0.5*(1-2*ux)*c_12_post + 0.5*c_22_post
    cDash_0_2_post  =          (1-ux*ux)*c_02_post -         2*ux*c_12_post -     c_22_post
    cDash_1_2_post  =  0.5*(ux*(1 + ux))*c_02_post + 0.5*(1+2*ux)*c_12_post + 0.5*c_22_post

    # post collision distributions
    f__n1_n1_post = -0.5*(uy*(1 - uy))*cDash_n1_0_post - 0.5*(1-2*uy)*cDash_n1_1_post + 0.5*cDash_n1_2_post
    f__n1_0_post  =        (1 - uy*uy)*cDash_n1_0_post -         2*uy*cDash_n1_1_post -     cDash_n1_2_post
    f__n1_1_post  = -0.5*(uy*(1 + uy))*cDash_n1_0_post + 0.5*(1+2*uy)*cDash_n1_1_post + 0.5*cDash_n1_2_post

    f__0_n1_post = -0.5*(uy*(1 - uy))*cDash_0_0_post - 0.5*(1-2*uy)*cDash_0_1_post + 0.5*cDash_0_2_post
    f__0_0_post  =        (1 - uy*uy)*cDash_0_0_post -         2*uy*cDash_0_1_post -     cDash_0_2_post
    f__0_1_post  = -0.5*(uy*(1 + uy))*cDash_0_0_post + 0.5*(1+2*uy)*cDash_0_1_post + 0.5*cDash_0_2_post

    f__1_n1_post = -0.5*(uy*(1 - uy))*cDash_1_0_post - 0.5*(1-2*uy)*cDash_1_1_post + 0.5*cDash_1_2_post
    f__1_0_post  =        (1 - uy*uy)*cDash_1_0_post -         2*uy*cDash_1_1_post -     cDash_1_2_post
    f__1_1_post  = -0.5*(uy*(1 + uy))*cDash_1_0_post + 0.5*(1+2*uy)*cDash_1_1_post + 0.5*cDash_1_2_post

    return array([f__0_0_post, f__1_0_post, f__0_1_post, f__n1_0_post, f__0_n1_post,   f__1_1_post, f__n1_1_post, f__n1_n1_post, f__1_n1_post])
