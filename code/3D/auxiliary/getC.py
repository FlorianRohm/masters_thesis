from numpy import zeros

def getC():
    c = zeros((3, 3, 3, 3), dtype='int64')
    # dimensions of c are
    # the three indices of D3Q27, ijk
    # the velocity vector for the index
    for i in range(3):
        for j in range(3):
            for k in range(3):
                c[i, j, k, 0] = (i-1)
                c[i, j, k, 1] = (j-1)
                c[i, j, k, 2] = (k-1)
    return c
