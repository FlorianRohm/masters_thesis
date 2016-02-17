#!/usr/bin/python
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License, either
# version 3 of the License, or (at your option) any later version.
#
#
# 2D flow in a periodic box
#
# D2Q9 Stencil with enumeration
#
# 6   2   5
#   \ | /
# 3 - 0 - 1
#   / | \
# 7   4   8
#

from numpy import *
from matplotlib import pyplot
from auxiliary.collide import BGKCollide, cumulantCollide, cumulantCollideAll, centralMomentSRT
from auxiliary.stream import stream
from auxiliary.LBMHelpers import getMacroValues,equilibrium
from auxiliary.visualize import visualize
from auxiliary.consistency import checkTransformation

import os

###### Plot settings ############################################################

plotEveryN    = 1         # draw every plotEveryN'th cycle
skipFirstN    = 0       # do not process the first skipFirstN cycles
savePlot      = True      # save velocity norm and x velocity plot
liveUpdate    = True      # show the process of the simulation (slow)
saveVTK       = False       # save the vtk files
prefix        = 'cylinder'      # naming prefix for saved files
outputFolder  = './out'    # folder to save the output to
workingFolder = os.getcwd()


###### Flow definition #########################################################
maxIterations = 230  # Total number of time iterations.
Re            = 100.0   # Reynolds number.re

# Number of Cells
ny = 100
nx = 100

# Highest index in each direction
nxl = nx-1
nyl = ny-1

# populations
q  = 9

# Velocity in lattice units.
uLB  = 0.04
uLBAverage = uLB
nulb = uLBAverage*nx/Re

# Relaxation parameter
omega = 1.0 / (3.*nulb+0.5)

print "Omega: {}".format(omega)

###### Plot preparations ############################################################

# quick and dirty way to create output directory
if not os.path.isdir(outputFolder):
    try:
        os.makedirs(outputFolder)
    except OSError:
        pass

# Define the Grid for vtk output
gridX = arange(0, nx, dtype='float64')
gridY = arange(0, ny, dtype='float64')
gridZ = arange(0, 1, dtype='float64')
grid  = gridX, gridY, gridZ

# Set velocity to 0 on z-axis
velocityZ = zeros((nx, ny, 1))

# define axis for velocity plot
axisYPlot = arange(0, ny, dtype='float64')
axisYNorm = axisYPlot/(ny)


# constant start velocity uLB in x direction for the lower half
delimiter = fromfunction(lambda x, y: (y*2 + sin(x/10)*5 < nxl),  (nx, ny))
vel = fromfunction(lambda d, x, y: (1-d)*uLB*2,  (2, nx, ny))
vel[0, delimiter] = 0;

feq   = equilibrium(1.0, vel)
fin   = feq.copy()

# interactive mode (execute code while showing figures)
if ( liveUpdate ):
    pyplot.ion()
    fig, ax = pyplot.subplots(1)

os.chdir(outputFolder)

plottingData = (plotEveryN, skipFirstN, liveUpdate, saveVTK, savePlot, ax, fig, grid, prefix)

checkTransformation(fin)

###### Main time loop ##########################################################
for time in range(maxIterations):

    # Calculate macroscopic density and velocity
    (rho, u) = getMacroValues(fin)

    predensity = sum(rho)
    preVelX = sum(u[0,:,:])
    preVelY = sum(u[1,:,:])


    feq = equilibrium(rho, u)

    # Collision step.
    #fpost = BGKCollide(fin, feq, omega)
    #fpost = cumulantCollide(fin, rho, u, omega)
    #fpost = cumulantCollideAll(fin, rho, u, omega, omega, omega, omega)
    fpost = centralMomentSRT(fin, feq, u, omega)

    (rho, u) = getMacroValues(fpost)

    deltaRho = sum(rho) - predensity
    deltaVelX = sum(u[0,:,:]) - preVelX
    deltaVelY = sum(u[0,:,:]) - preVelY
    maxRho = amax(rho)
    minRho = amin(rho)
    print minRho

    # Streaming step
    fin = stream(fpost)

    visualize(u, rho, time, *plottingData)

os.chdir(workingFolder)
