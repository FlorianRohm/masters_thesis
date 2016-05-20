from numpy import *

def stream(f):
    post = f.copy()
    post[0, :, :] = f[0, :, :]

    post[1, :, :] = roll(f[1, :, :],  1,  axis=0)
    post[2, :, :] = roll(f[2, :, :],  1,  axis=1)
    post[3, :, :] = roll(f[3, :, :], -1,  axis=0)
    post[4, :, :] = roll(f[4, :, :], -1,  axis=1)

    post[5, :, :] = roll(roll(f[5, :, :],  1,  axis=0),  1,  axis=1)
    post[6, :, :] = roll(roll(f[6, :, :], -1,  axis=0),  1,  axis=1)
    post[7, :, :] = roll(roll(f[7, :, :], -1,  axis=0), -1,  axis=1)
    post[8, :, :] = roll(roll(f[8, :, :],  1,  axis=0), -1,  axis=1)

    return post

def stream2(f):
    post = f.copy()
    (n, nx, ny) = f.shape
    nxl = nx-1
    nyl = ny-1

    post[0, :, :] = f[0, :, :]

    post[1, 1:nxl,   :]     = f[1, 0:nxl-1,  :]
    post[2,   :,   0:nyl-1] = f[2,   :,    1:nyl]
    post[3, 0:nxl-1, :]     = f[3, 1:nxl,    :]
    post[4,   :,   1:nyl]   = f[4,   :,    0:nyl-1]

    post[5, 1:nxl,   0:nyl-1] = f[5, 0:nxl-1, 1:nyl]
    post[6, 0:nxl-1, 0:nyl-1] = f[6, 1:nxl,   1:nyl]
    post[7, 0:nxl-1, 1:nyl]   = f[7, 1:nxl,   0:nyl-1]
    post[8, 1:nxl,   1:nyl]   = f[8, 0:nxl-1, 0:nyl-1]

    post[1, 0,   :] = f[1, nxl,  :]
    post[2, :, nyl] = f[2,  :,   0]
    post[3, nxl, :] = f[3, 0,    :]
    post[4, :,   0] = f[4,  :, nyl]

    post[5, 0,   nyl] = f[5, nxl, 0]
    post[6, nxl, nyl] = f[6, 0,   0]
    post[7, nxl, 0]   = f[7, 0,   nyl]
    post[8, 0,   0]   = f[8, nxl, nyl]

    return post
