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
# 5   2   8
#   \ | /
# 3 - 0 - 6
#   / | \
# 4   1   7
#

from numpy import *
from matplotlib import cm, pyplot
from auxiliary.VTKWrapper import saveToVTK
from auxiliary.stream import stream
from auxiliary.collide import BGKCollide, cumulantCollide, cumulantCollideAllInOne
from auxiliary.LBMHelpers import clamp, getMacroValues, sumPopulations, equilibrium, noslip, iLeft, iCentV, iRight, iTop, iCentH, iBot
from auxiliary.ghiaResults import *
from auxiliary.transformations.momentsFromDistributions import momentsFromDistributions
import os


###### Plot settings ############################################################

plotEveryN    = 100         # draw every plotEveryN'th cycle
skipFirstN    = 0       # do not process the first skipFirstN cycles
savePlot      = False      # save velocity norm and x velocity plot
liveUpdate    = True     # show the process of the simulation (slow)
saveVTK       = False       # save the vtk files
prefix        = 'ldc_20k'      # naming prefix for saved files
outputFolder  = './out'    # folder to save the output to
workingFolder = os.getcwd()


###### Flow definition #########################################################
maxIterations = 200000  # Total number of time iterations.
Re            = 200.0   # Reynolds number.re

# Number of Cells
nx = 100
ny = 100

# Highest index in each direction
nxl = nx-1
nyl = ny-1

# populations
q  = 9

# Velocity in lattice units.
uLB  = 0.08
nulb = uLB*ny/Re

# Relaxation parameter
omega = 1.0 / (3.*nulb+0.5)


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
axisYPlot = arange(ny, 0, -1, dtype='float64')
axisYNorm = axisYPlot/(ny)

# import results of Ghia et al.
yghia = getAxis()
xghia = getRe100()


###### Setup ##################################################################

# set up the walls
leftWall    = fromfunction(lambda   x, y: x == 0,  (nx, ny))
rightWall   = fromfunction(lambda   x, y: x == nxl,  (nx, ny))
bottomWall  = fromfunction(lambda   x, y: y == nyl,  (nx, ny))

solidDomain = logical_or(logical_or(leftWall, rightWall), bottomWall)
fluidDomain = invert(solidDomain)

# initial velocity profile
#  < -
#   -
#  - >
vel                = zeros((2, nx, ny))
vel[0, :,  0 ]     =  uLB

# initial particle distributions
feq   = equilibrium(1.0, vel.reshape((2,nx*ny))).reshape((9,nx,ny))
fin   = feq.copy()
fpost = feq.copy()  # post collision distributions

# interactive mode (execute code while showing figures)
if ( liveUpdate | savePlot ):
    pyplot.ion()
    fig, ax = pyplot.subplots(1)

os.chdir(outputFolder)


###### Main time loop ##########################################################
for time in range(maxIterations):

    (rho, u) = getMacroValues(fin)

    # Collision step.
    #fpost = BGKCollide(fin, feq, omega)
    fpost[:,fluidDomain] = cumulantCollideAllInOne(fin[:,fluidDomain], rho[fluidDomain], u[:,fluidDomain], omega )

    # Streaming step
    fin = stream(fpost)

    # Left wall: compute density from known populations
    u[:, :, 0] = vel[:, :, 0]
    rho[:, 0] = sumPopulations(fin[iCentV, 0, :])+2.*sumPopulations(fin[iTop, 0, :])
    # Calculate macroscopic density and velocity

    feq[:,:,0] = equilibrium(rho[:,0], u[:,:,0])

    # complete the left wall treatement wrt Yu 2002
    fin[iBot, 0, :] = - feq[iTop, 0, :] + (feq[iBot, 0, :] + fin[iTop, 0, :])
    # fin[iBot, 0, :] = fin[iTop, 0, :] + 6 * dot(c, u.transpose(1, 0, 2)) c[iBot][0] * rho*t[iBot]

    # bounce back distributions at walls
    finc = fin[:, solidDomain].copy()
    #   print finc.shape
    for i in range(9):
        fin[i, solidDomain] = finc[noslip[i], :]

    # Visualization
    if ( (time % plotEveryN == 0) & (liveUpdate | saveVTK | savePlot) & (time > skipFirstN) ):
        if ( liveUpdate | savePlot ):
            ax.clear()
            ax.imshow(sqrt(u[0]**2+u[1]**2).transpose(),  cmap=cm.jet, vmin=0., vmax=0.01)
            ax.set_title('velocity norm')

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

os.chdir(workingFolder)
