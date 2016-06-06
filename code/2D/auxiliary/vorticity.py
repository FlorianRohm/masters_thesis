from numpy import *

def vorticity(u):
    ux = u[0,:,:]
    uy = u[1,:,:]

    duxdy = (roll(ux,-1,axis=1) - roll(ux,1,axis=1))/2.
    duydx = (roll(uy,-1,axis=0) - roll(uy,1,axis=0))/2.

    return duydx - duxdy
