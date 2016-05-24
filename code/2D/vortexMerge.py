#!/usr/bin/python
# Copyright (C) 2013 FlowKit Ltd, Lausanne, Switzerland
# E-mail contact: contact@flowkit.com
#
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License, either
# version 3 of the License, or (at your option) any later version.
#
# Rewrite, extension and bugfixes by Florian Rohm, 2015

#
# 2D flow around a cylinder
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
from matplotlib import cm, pyplot
from auxiliary.VTKWrapper import saveToVTK
from auxiliary.collide import BGKCollide, cumulantCollideAllInOne
from auxiliary.vortex import rankinVortex
from auxiliary.stream import stream
from auxiliary.LBMHelpers import clamp, getMacroValues, sumPopulations, equilibrium, noslip, iLeft, iCentV, iRight, iTop, iCentH, iBot
from auxiliary.boundaryConditions import YuLeft
from auxiliary.obstacle import obstacleAttack, drag, lift
from random import random
import sys, getopt
import os
import datetime

###### Parsing ############################################################

Re = 100
size = 100
collisionFunction = BGKCollide
collStr = "srt"
try:
  opts, args = getopt.getopt(sys.argv[1:],"hcbr:s:",)
except getopt.GetoptError:
    print "parse error"
    print 'test.py -i <inputfile> -o <outputfile> -f <sampling frequency in diameters of sphere>'
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print 'test.py -i <inputfile> -o <outputfile> -f <sampling frequency>'
        sys.exit()
    elif opt in ("-r"):
        Re = int(float(arg))
    elif opt in ("-s"):
        size = int(float(arg))
    elif opt in ("-c"):
        collisionFunction = cumulantCollideAllInOne
        collStr = "cumulant"
    elif opt in ("-b"):
        collisionFunction = BGKCollide
        collStr = "srt"

plotEveryN = 150
print 'Begin of calculation with {0} collision with Re={1} and size {2}'.format(collStr,Re,size)
startTime = datetime.datetime.now()
factor = size/10.


###### Plot settings ############################################################

skipFirstN  = 0       # initial conditions already quite interesting
savePlot      = False      # save velocity norm and x velocity plot
liveUpdate    = True      # show the process of the simulation (slow)
saveVTK       = False       # write out drag and lift
prefix        = 'vortexMerge_{0}_Re{1}_size{2}'.format(collStr, Re, size)      # naming prefix for saved files
outputFolder  = './out'    # folder to save the outputFile to
workingFolder = os.getcwd()

###### Flow definition #########################################################
maxIterations = size*200  # Total number of time iterations.
maxIterations = 1000

# Number of Cells
ny = size
nx = size

# Highest index in each direction
nxl = nx-1
nyl = ny-1

# populations
q  = 9

# Velocity in lattice units.
uLB  = 0.04
uLBAverage = 2./3.*uLB # according to schaefer turek 2D-2
nulb = uLBAverage*size/Re

# Relaxation parameter
omega = 1.0 / (3.*nulb+0.5)

###### Plot preparations ############################################################

# quick and dirty way to create outputFile directory
if not os.path.isdir(outputFolder):
    try:
        os.makedirs(outputFolder)
    except OSError:
        pass
os.chdir(outputFolder)

# Define the Grid for vtk output
gridX = arange(0, nx, dtype='float64')
gridY = arange(0, ny, dtype='float64')
gridZ = arange(0, 1, dtype='float64')
grid  = gridX, gridY, gridZ

# Set velocity to 0 on z-axis
velocityZ = zeros((nx, ny, 1))

###### Setup ##################################################################

# velocity
radius1 = 5.
v1x = size/2
v1y = size/2 - 2*radius1
radius2 = 5.
v2x = size/2
v2y = size/2 +2*radius2

vel = zeros((2, nx, ny))
for x in range(nx):
    for y in range(ny):
        vel[:,x,y] = vel[:,x,y] + rankinVortex(x-v1x, y-v1y, radius1) + rankinVortex(x-v2x, y-v2y, radius2)

vel = vel*uLB
# initial particle distributions
feq   = equilibrium(1.0, vel)
fin   = feq.copy()
fpost = feq.copy()  # post collision distributions


# interactive mode (execute code while showing figures)
if ( liveUpdate | savePlot ):
    pyplot.ion()
    fig, ax = pyplot.subplots(1)
    ax.clear()
    ax.imshow(sqrt(vel[0]**2+vel[1]**2).transpose(),  cmap=cm.afmhot, vmin=0., vmax=0.05)


###### Main time loop ##########################################################
for time in range(maxIterations):
    # Calculate macroscopic density and velocity
    (rho, u) = getMacroValues(fin)

    # Collision step.
    fpost = collisionFunction(fin, rho, u, omega )

    # Streaming step
    fin = stream(fpost)

    # Visualization
    if ( (time % plotEveryN == 0) & (liveUpdate | saveVTK | savePlot) & (time > skipFirstN) ):
        if ( liveUpdate | savePlot ):
            ax.clear()
            ax.imshow(sqrt(u[0]**2+u[1]**2).transpose(),  cmap=cm.afmhot, vmin=0., vmax=0.05)

        if ( liveUpdate ):
            pyplot.draw()
        if ( saveVTK ):
            # convert velocity and density to 3d arrays
            printVel = reshape(u, (2, nx, ny, 1))
            printRho = reshape(rho, (nx, ny, 1))

            velocity = (printVel[0, :, :, :], printVel[1, :, :, :], velocityZ)
            saveNumber = str(time/plotEveryN).zfill(4)

            saveToVTK(velocity, printRho, prefix, saveNumber, grid)
        if ( savePlot ):
            pyplot.savefig(prefix + "." + str(time/plotEveryN).zfill(4) + ".png")

endTime = datetime.datetime.now()
deltaTime = endTime - startTime
timeFile = open("{0}_time".format(prefix),"w")
timeFile.write("{0}".format(deltaTime.total_seconds()));
timeFile.write("\n");
timeFile.close()
os.chdir(workingFolder)
print 'End of calculation with {0} collision with Re={1} and size {2}.\n Elapsed time: {3}'.format(collStr,Re,size,deltaTime.total_seconds())
