from numpy import *
from equilibrium import equilibrium

def BGKCollide(fin, rho, u, omega):
    feq = equilibrium(rho, u)
    return fin - omega * (fin - feq)
