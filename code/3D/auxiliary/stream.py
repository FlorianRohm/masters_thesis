from numpy import roll
from getC import getC

c = getC()

def stream(f):
    for i in range(3):
        for j in range(3):
            for k in range(3):
                roll(roll(roll( f[i,j,k], c[i,j,k,0], axis=0), c[i,j,k,1], axis=1), c[i,j,k,2], axis=2)
    return f
