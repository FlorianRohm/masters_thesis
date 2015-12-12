from numpy import *


###### Lattice Constants #######################################################
q = 9

# lattice velocities
c = array([[0, 0], [1, 0], [0, 1], [-1, 0], [0, -1], [1, 1], [-1, 1], [-1, -1], [1, -1]])

# Lattice weights
t      = 1./36. * ones(q)
t[1:5] = 1./9.
t[0]   = 4./9.

# index array for noslip reflection
# inverts the d2q9 stencil for use in bounceback scenarios
noslip = [0, 3, 4, 1, 2, 7, 8, 5, 6]

# index arrays for different sides of the d2q9 stencil
iLeft   = arange(q)[asarray([ci[0] <  0 for ci in c])]
iCentV  = arange(q)[asarray([ci[0] == 0 for ci in c])]
iRight  = arange(q)[asarray([ci[0] >  0 for ci in c])]
iTop    = arange(q)[asarray([ci[1] >  0 for ci in c])]
iCentH  = arange(q)[asarray([ci[1] == 0 for ci in c])]
iBot    = arange(q)[asarray([ci[1] <  0 for ci in c])]


###### Function Definitions ####################################################

# Helper function for density computation.
def sumPopulations(fin):
    return sum(fin, axis = 0)


# Equilibrium distribution function.
def equilibrium(rho, u):
    (d, nx, ny) = u.shape
    cu   = 3.0 * dot(c, u.transpose(1, 0, 2))
    usqr = 3./2.*(u[0, :, :]**2+u[1, :, :]**2)
    feq = zeros((q, nx, ny))
    for i in range(q):
        feq[i, :, :] = rho*t[i]*(1.+cu[i]+0.5*cu[i]**2-usqr)
    return feq


def clamp(val, minVal, maxVal):
    return maximum(minVal, minimum(val, maxVal))

def getMacroValues(f):
    # Calculate macroscopic density ...
    rho = sumPopulations(f)
    # ... and velocity
    u = dot(c.transpose(), f.transpose((1, 0, 2)))/rho
    return (rho, u)
