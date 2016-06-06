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
from random import random, triangular, seed
from numpy.fft import fft2
from auxiliary.vorticity import vorticity
seed(1503)

import sys, getopt
import os
import datetime

###### Parsing ############################################################

Re = 100
size = 128
collisionFunction = BGKCollide
collStr = "srt"
try:
  opts, args = getopt.getopt(sys.argv[1:],"hcbr:s:",)
except getopt.GetoptError:
    print "parse error"
    print 'test.py -s <size of grid> -r <reynolds number> -c <use cumulant collision> '
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print 'test.py -s <size of grid> -r <reynolds number> -c <use cumulant collision> '
        sys.exit()
    elif opt in ("-r"):
        Re = int(float(arg))
    elif opt in ("-s"):
        size = int(float(arg))
    elif opt in ("-c"):
        collisionFunction = cumulantCollideAllInOne
        collStr = "cumulant"

plotEveryN = 100
print 'Begin of calculation with {0} collision with Re={1} and size {2}'.format(collStr,Re,size)
startTime = datetime.datetime.now()
factor = size/10.


###### Plot settings ############################################################

skipFirstN  = 0       # initial conditions already quite interesting
savePlot      = True      # save velocity norm and x velocity plot
liveUpdate    = True      # show the process of the simulation (slow)
analysis      = False       # calculate and write enstrophy
saveVTK       = False       # write out drag and lift
prefix        = '{0}_Re{1}_size{2}'.format(collStr, Re, size)      # naming prefix for saved files
outputFolder  = './out/vortexMerge'    # folder to save the outputFile to
workingFolder = os.getcwd()

###### Flow definition #########################################################
maxIterations = size*4000  # Total number of time iterations.

# Number of Cells
ny = size
nx = size

# Highest index in each direction
nxl = nx-1
nyl = ny-1

# populations
q  = 9

# Velocity in lattice units.
uLB  = 0.01

diameter = 10
nulb = uLB*diameter/Re

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
c1x = size/2
c1y = size/2 - diameter
c2x = size/2
c2y = size/2 + diameter


# velocity
print "Building vortices"
vel = zeros((2, nx, ny))
for x in range(nx):
    for y in range(ny):
        vel[:,x,y] = vel[:,x,y] \
        + rankinVortex(x-c1x, y-c1y, diameter/2.) \
        + rankinVortex(x-c2x, y-c2y, diameter/2.) \

vel = vel*uLB
print "Building vortices complete"
avgVelX = mean(vel[0,:,:])
avgVelY = mean(vel[1,:,:])

# initial particle distributions
feq   = equilibrium(1.0, vel.reshape((2,nx*ny))).reshape((9,nx,ny))
fin   = feq.copy()
fpost = feq.copy()  # post collision distributions

fluidDomain = full((nx,ny), True, dtype=bool)
clipfft = 30
fftDomain = fluidDomain[nx/2-clipfft:nx/2+clipfft, ny/2-clipfft:ny/2+clipfft]

# interactive mode (execute code while showing figures)
if ( liveUpdate | savePlot ):
    pyplot.ion()
    fig, ax = pyplot.subplots(1)
    ax.clear()
    ax.imshow(sqrt(vel[0]**2+vel[1]**2).transpose(),  cmap=cm.afmhot, vmin=0., vmax=0.05)

if analysis:
    outputFile = open(prefix, 'w')

###### Main time loop ##########################################################
for time in range(maxIterations):
    # Calculate macroscopic density and velocity
    (rho, u) = getMacroValues(fin)

    # Collision step.
    fpost[:,fluidDomain] = collisionFunction(fin[:,fluidDomain], rho[fluidDomain], u[:,fluidDomain], omega )

    # Streaming step
    fin = stream(fpost)

    # Visualization
    if ( (time % plotEveryN == 0) & (liveUpdate | saveVTK | savePlot | analysis) & (time > skipFirstN) ):
        vort = vorticity(u);
        meanEnstrophy = mean(vort**2)/2.
        turbulentKineticEnergy = mean((u[0]-avgVelX)**2 + (u[1]-avgVelY)**2)/2.

        if (analysis):
            outputFile.write("{0},{1}\n".format(meanEnstrophy, turbulentKineticEnergy))

        if ( liveUpdate | savePlot ):
            ax.clear()
            #vel2 = u[0]**2+u[1]**2
            #vel2 = vel2.transpose();
            #fftVel = fft2(vel2)
            #fftVel = roll(roll(abs(fftVel),nx/2, axis=0), ny/2, axis=1)
            #fftVel = fftVel[fftDomain].reshape((2*clipfft,2*clipfft))

            ax.imshow(vort,  cmap=cm.gray, vmin=-0.005, vmax=0.005 )
            #ax.imshow(sqrt(u[0]**2+u[1]**2).transpose(),  cmap=cm.afmhot, vmin=0., vmax=0.05)

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
if analysis:
    timeFile = open("{0}_time".format(prefix),"w")
    timeFile.write("{0},{1}".format(deltaTime.total_seconds(), sum(fluidDomain)));
    timeFile.write("\n");
    timeFile.close()
    outputFile.close()
os.chdir(workingFolder)
print 'End of calculation with {0} collision with Re={1} and size {2}.\n Elapsed time: {3}'.format(collStr,Re,size,deltaTime.total_seconds())
