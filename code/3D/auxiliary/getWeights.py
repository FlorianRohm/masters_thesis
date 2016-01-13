from numpy import ones

def getWeights():
    t = ones((3,3,3))/54. # edge centers
    t[1,1,1] = 8./27. # center
    t[0,1,1] = t[1,0,1] = t[1,1,0] = \
    t[2,1,1] = t[1,2,1] = t[1,1,2] = 2./27.  # face centers

    t[0,0,0] = t[0,0,2] = t[0,2,0] = t[0,2,2] = \
    t[2,0,0] = t[2,0,2] = t[2,2,0] = t[2,2,2] = 1./216. # corners

    return t
