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
    (d, nx) = u.shape
    cu   = 3.0 * dot(c, u[:,:])
    usqr = 3./2.*(u[0, :]**2+u[1, :]**2)
    feq = zeros((q, nx))
    rhot1 = t[1]*rho
    rhot5 = t[5]*rho
    feq[1, :] = rhot1*(1.+cu[1]+0.5*cu[1]**2-usqr)
    feq[2, :] = rhot1*(1.+cu[2]+0.5*cu[2]**2-usqr)
    feq[3, :] = rhot1*(1.+cu[3]+0.5*cu[3]**2-usqr)
    feq[4, :] = rhot1*(1.+cu[4]+0.5*cu[4]**2-usqr)
    feq[5, :] = rhot5*(1.+cu[5]+0.5*cu[5]**2-usqr)
    feq[6, :] = rhot5*(1.+cu[6]+0.5*cu[6]**2-usqr)
    feq[7, :] = rhot5*(1.+cu[7]+0.5*cu[7]**2-usqr)
    feq[8, :] = rhot5*(1.+cu[8]+0.5*cu[8]**2-usqr)
    feq[0, :] = rho*t[0]*(1.-usqr)
    return feq


def clamp(val, minVal, maxVal):
    return maximum(minVal, minimum(val, maxVal))

def getMacroValues(f):
    # Calculate macroscopic density ...
    rho = sumPopulations(f)
    # ... and velocity
    u = dot(c.transpose(), f.transpose((1, 0, 2)))/rho
    return (rho, u)
