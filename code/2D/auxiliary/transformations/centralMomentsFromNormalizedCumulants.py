
def centralMomentsFromNormalizedCumulants (rho, K_00, K_10, K_01, K_11, K_20, K_02, K_21, K_12, K_22):
    # Transformation to central moments
    c_22 = K_22 + 2*K_11*K_11/rho + K_20*K_02/rho

    return (rho, 0, 0, K_11, K_20, K_02, K_21, K_12, c_22)
