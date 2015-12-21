from numpy import *

def BGKCollide(fin, feq, omega):
    return fin - omega * (fin - feq)


def cumulantCollide(fin, omega, u):
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
    K_22 = c_22 - 2*c_11*c_11 - c_20*c_02

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
    c_22_post = K_22_post + 2*c_11_post*c_11_post + c_20_post*c_02_post


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
