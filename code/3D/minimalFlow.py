#!/usr/bin/python

# 3D Flow in a periodic box
# using the D3Q27 stencil
#
# z    y
# |   /
# |  /
# | /
#    ------- x


from numpy import *

from auxiliary.equilibrium import equilibrium
from auxiliary.getMacroValues import getMacroValues
from auxiliary.collide import BGKCollide, cumulantCollide
from auxiliary.stream import stream
from auxiliary.visualize import visualize

import os

###### settings ############################################################

saveEveryN    = 5
skipFirstN    = 0
prefix        = 'minimal'      # naming prefix for saved files
outputFolder  = './out'    # folder to save the output to
workingFolder = os.getcwd()

###### Flow definitions ####################################################

maxIterations = 5
Re            = 100. # Reynolds Number

# Number of Cells
nx = 100
ny = 3
nz = 100

# highest index per direction
nxl = nx -1
nyl = ny -1
nzl = nz -1

# Velocity in lattice units
uLB = 0.04
nulb = uLB*nx/Re

omega = 1. / (3*nulb + 0.5)

print "Omega: {}".format(omega)

###### Plot preparations ###################################################

# quick and dirty way to create output directory
if not os.path.isdir(outputFolder):
    try:
        os.makedirs(outputFolder)
    except OSError:
        pass

# Define the Grid for vtk output
gridX = arange(0, nx, dtype='float64')
gridY = arange(0, ny, dtype='float64')
gridZ = arange(0, nz, dtype='float64')
grid  = gridX, gridY, gridZ

###### Flow preparations ####################################################

# constant start velocity uLB in x direction for the lower half
delimiter = fromfunction(lambda x, y, z: (z*2 + sin(x/10)*3 + sin(y/10)*3 < nzl),  (nx, ny, nz))
vel = fromfunction(lambda d, x, y, z: (1-d)*(2-d)*uLB*2,  (3, nx, ny, nz))
vel[0, delimiter] = 0;

fin = equilibrium(1.0, vel)

os.chdir(outputFolder)

savingOptions = (saveEveryN, skipFirstN, grid, prefix)
###### Main time loop ##########################################################
for time in range(maxIterations):

    (rho, u) = getMacroValues(fin)
    #minRho = amin(rho)
    #print minRho

    #fpost = BGKCollide(fin, rho, u, omega)
    fpost = cumulantCollide(fin, rho, u, omega)

    fin = stream(fpost)

    visualize(rho, u, time, *savingOptions)

os.chdir(workingFolder)
