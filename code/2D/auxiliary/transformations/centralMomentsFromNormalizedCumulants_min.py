
def centralMomentsFromNormalizedCumulants_min (rho, K_11, K_20, K_02):
    # Transformation to central moments
    c_22 = 2*K_11*K_11/rho + K_20*K_02/rho

    return (rho, K_11, K_20, K_02, c_22)
