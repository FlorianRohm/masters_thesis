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
from auxiliary.collide import BGKCollide, cumulantCollide
from auxiliary.stream import stream
from auxiliary.LBMHelpers import clamp, getMacroValues, sumPopulations, equilibrium, noslip, iLeft, iCentV, iRight, iTop, iCentH, iBot
from auxiliary.boundaryConditions import YuLeft
from auxiliary.obstacle import obstacleAttack, drag, lift

import os


###### Plot settings ############################################################

plotEveryN    = 100         # draw every plotEveryN'th cycle
skipFirstN    = 0       # do not process the first skipFirstN cycles
savePlot      = True      # save velocity norm and x velocity plot
liveUpdate    = False      # show the process of the simulation (slow)
saveVTK       = False       # save the vtk files
prefix        = 'cyl_Re5k'      # naming prefix for saved files
outputFolder  = './out'    # folder to save the output to
workingFolder = os.getcwd()


###### Flow definition #########################################################
maxIterations = 200000  # Total number of time iterations.
Re            = 5000.0   # Reynolds number.re

# Number of Cells
ny = 100
nx = 300

# Highest index in each direction
nxl = nx-1
nyl = ny-1

# populations
q  = 9

# Coordinates of the cylinder.
cx = ny*0.0909/0.18636
cy = ny*0.512
r  = ny/44./0.18636
diameter = 2*r

# Velocity in lattice units.
uLB  = 0.04
uLBAverage = 2./3.*uLB
nulb = uLBAverage*diameter/Re

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
axisYPlot = arange(0, ny, dtype='float64')
axisYNorm = axisYPlot/(ny)

###### Setup ##################################################################

# cylindrical obstacle
obstacle       = fromfunction(lambda    x, y: (x-cx)**2+(y-cy)**2 < r**2,  (nx, ny))
obstacleBounds, dragBoundStencil, liftBoundStencil, completeBoundStencil = obstacleAttack(obstacle)

boundary = fromfunction(lambda x, y: logical_or((y == 0), (y == ny)), (nx, ny))

# velocity inlet with small perturbation
# vel = fromfunction(lambda d, x, y: (1-d)*uLB*(1.0+1e-2*sin(y/nyl*2*pi)),  (2, nx, ny))

# velocity inlet for schaefer turek
vel = fromfunction(lambda d, x, y: (1-d)*4*uLB*y*(nyl-y)/(nyl**2),  (2, nx, ny))

# initial particle distributions
feq   = equilibrium(1.0, vel)
# feq   = equilibrium(1.0, zeros((2, nx, ny)))
fin   = feq.copy()
fpost = feq.copy()  # post collision distributions


# interactive mode (execute code while showing figures)
if ( liveUpdate | savePlot ):
    pyplot.ion()
    fig, ax = pyplot.subplots(1)

os.chdir(outputFolder)


###### Main time loop ##########################################################
for time in range(maxIterations):
    # bounce back distributions at obstacle
    for i in range(q):
        fin[i, obstacle] = fin[noslip[i], obstacle]
    # and walls
    for i in range(q):
        fin[i, boundary] = fin[noslip[i], boundary]

    # Right Wall: Produce zero pressure gradient for the outflow
    fin[iLeft, -1, :] = fin[iLeft, -2, :]

    # Calculate macroscopic density and velocity
    (rho, u) = getMacroValues(fin)

    # Left wall: compute density from known populations.
    u[:, 0, :] = vel[:, 0, :]
    rho[0, :] = 1./(1.-u[0, 0, :]) * (sumPopulations(fin[iCentV, 0, :])+2.*sumPopulations(fin[iLeft, 0, :]))

    feq = equilibrium(rho, u)

    # complete the left wall treatement wrt Yu 2002
    fin[iRight, 0, :] = feq[iLeft, 0, :] + (feq[iRight, 0, :] - fin[iLeft, 0, :])

    # Collision step.
    #fpost = BGKCollide(fin, feq, omega)
    fpost = cumulantCollide(fin, rho, u, omega)

    # Streaming step
    fin = stream(fpost)

    # Visualization
    if ( (time % plotEveryN == 0) & (liveUpdate | saveVTK | savePlot) & (time > skipFirstN) ):
        # Here, distributions are streamed into the obstacle -> compute drag and lift
        scaling = 2/(sumPopulations(fin)*uLBAverage*uLBAverage*diameter)
        scaledFin = fin.copy()
        # just scale the populations which are used
        for i in range(9):
            scaledFin[i, completeBoundStencil] = scaledFin[i, completeBoundStencil]*scaling[completeBoundStencil]


        dragCoeff = drag(scaledFin, obstacleBounds)
        liftCoeff = lift(scaledFin, obstacleBounds)

        if ( liveUpdate | savePlot ):
            ax.clear()
            ax.imshow(sqrt(u[0]**2+u[1]**2).transpose(),  cmap=cm.afmhot, vmin=0., vmax=0.1)
            ax.set_title('velocity norm')

            textstr = '$\mathrm{drag}= %.4f$\n$\mathrm{lift}= %.4f$' % (dragCoeff, liftCoeff)
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            ax.text(0.95, 0.95,
                    textstr,
                    horizontalalignment='right', verticalalignment='top',
                    transform=ax.transAxes, fontsize=14,
                    bbox=props)

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
