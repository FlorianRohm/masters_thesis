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
from auxiliary.boundary import cumulantBoundary
from auxiliary.stream import stream
from auxiliary.collide import BGKCollide, cumulantCollideAllInOne
from auxiliary.LBMHelpers import clamp, getMacroValues, sumPopulations, equilibrium, noslip, iLeft, iCentV, iRight, iTop, iCentH, iBot
from auxiliary.boundaryConditions import YuLeft
from auxiliary.obstacle import obstacleAttack, drag, lift
import sys, getopt
import os
import datetime

###### Parsing ############################################################

Re = 100
size = 10
frequency = 1         # draw every frequency'th cycle
collisionFunction = BGKCollide
collStr = "srt"

try:
  opts, args = getopt.getopt(sys.argv[1:],"hcr:s:f:",)
except getopt.GetoptError:
    print "parse error"
    print 'test.py -s <size of sphere> -r <reynolds number> -c <use cumulant collision> -f <sampling frequency in diameters of sphere>'
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print 'test.py -s <size of sphere> -r <reynolds number> -c <use cumulant collision> -f <sampling frequency in diameters of sphere>'
        sys.exit()
    elif opt in ("-r"):
        Re = int(float(arg))
    elif opt in ("-s"):
        size = int(float(arg))
    elif opt in ("-f"):
        frequency = int(float(arg))
    elif opt in ("-c"):
        collisionFunction = cumulantCollideAllInOne
        collStr = "cumulant"

plotEveryN = frequency*size
print 'Begin of calculation with {0} collision with Re={1} and size {2}'.format(collStr,Re,size)
startTime = datetime.datetime.now()

factor = size/10.


###### Plot settings ############################################################

savePlot      = False      # save velocity norm and x velocity plot
liveUpdate    = True      # show the process of the simulation (slow)
analysis      = False       # write out drag and lift
prefix        = '{0}_Re{1}_size{2}'.format(collStr, Re, size)      # naming prefix for saved files
outputFolder  = './out/schaeferTurekMB'    # folder to save the outputFile to
workingFolder = os.getcwd()

###### Flow definition #########################################################
maxIterations = size*3000  # Total number of time iterations.
plotEveryN = 50

# Number of Cells
ny = int(round(41*factor)) + 2 # for boundary
nx = int(round(220*factor)) + 2 # for boundary
nx=nx/2
skipFirstN = nx*3       # do not process the first skipFirstN cycles, 3 times nx for
skipFirstN = 0

# Highest index in each direction
nxl = nx-1
nyl = ny-1

# populations
q  = 9

# Coordinates of the cylinder.
cx = int(round(20*factor)) + 1 # +1 for boundary
cy = int(round(21*factor)) + 1
r = size/2

# Velocity in lattice units.
uLB  = 0.04
uLBAverage = 2./3.*uLB # according to schaefer turek 2D-2
nulb = uLBAverage*size/Re

preComputeFactorForScaling = 2/(uLBAverage*uLBAverage*size)

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
if analysis:
    outputFile = open(prefix, 'w')
    outputFile.write("timestep,drag,lift\n")

###### Setup ##################################################################

# cylindrical obstacle
obstacle       = fromfunction(lambda    x, y: (x-cx)**2+(y-cy)**2 < r**2,  (nx, ny))
obstacleBounds, dragBoundStencil, liftBoundStencil, completeBoundStencil = obstacleAttack(obstacle)
nrOfDragBoundary = sum(obstacleBounds[1]) + sum(obstacleBounds[5]) + sum(obstacleBounds[8]) + \
    sum(obstacleBounds[3]) + sum(obstacleBounds[6]) + sum(obstacleBounds[7])
nrOfLiftBoundary = sum(obstacleBounds[2]) + sum(obstacleBounds[6]) + sum(obstacleBounds[5]) + \
    sum(obstacleBounds[4]) + sum(obstacleBounds[7]) + sum(obstacleBounds[8])

noSlipBoundary = fromfunction(lambda x, y: logical_or((y == 0), (y == ny)), (nx, ny))


# velocity inlet for schaefer turek
velIn = fromfunction(lambda d, x, y: (1-d)*4*uLB*y*(nyl-y)/(nyl**2),  (2, nx, ny-2))
vel = zeros((2,nx,ny));
vel[:,:,1:ny-1] = velIn

# initial particle distributions
feq   = equilibrium(1.0, vel)
fin   = feq.copy()
fpost = feq.copy()  # post collision distributions


# interactive mode (execute code while showing figures)
if ( liveUpdate | savePlot ):
    pyplot.ion()
    fig, ax = pyplot.subplots(1)


###### Main time loop ##########################################################
for time in range(maxIterations):
    # bounce back distributions at obstacle
    for i in range(q):
        fin[i, obstacle] = fin[noslip[i], obstacle]
    # and walls
    for i in range(q):
        fin[i, noSlipBoundary] = fin[noslip[i], noSlipBoundary]

    # Calculate macroscopic density and velocity
    (rho, u) = getMacroValues(fin)

    #right wall: pressure equal to one
    fin[:,nxl,:] = cumulantBoundary(1, u[:,nxl-1,:])

    #left wall:
    fin[:,0,:] = cumulantBoundary(rho[0,:], vel[:,0,:])

    # Collision step.
    fpost[:,1:nx-1,1:ny-1] = collisionFunction(fin[:,1:nx-1,1:ny-1], rho[1:nx-1,1:ny-1], u[:,1:nx-1,1:ny-1], omega )

    # Streaming step
    fin = stream(fpost)

    # Visualization
    if ( (time % plotEveryN == 0) & (analysis | liveUpdate  | savePlot) & (time > skipFirstN) ):
        # Here, distributions are streamed into the obstacle -> compute drag and lift
        scaling = preComputeFactorForScaling/sumPopulations(fin[:,completeBoundStencil])
        scaledFin = fin.copy()
        # just scale the populations which are used
        for i in range(9):
            scaledFin[i, completeBoundStencil] = scaledFin[i, completeBoundStencil]*scaling

        dragCoeff = drag(scaledFin, obstacleBounds)
        liftCoeff = lift(scaledFin, obstacleBounds)

        if ( liveUpdate | savePlot ):
            ax.clear()
            velocityMag =sqrt(u[0]**2+u[1]**2)
            velocityMag[obstacle] = NAN
            ax.imshow(velocityMag.transpose(),  cmap=cm.afmhot, vmin=0., vmax=0.1)
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

        if ( savePlot ):
            pyplot.savefig(prefix + "." + str(time/plotEveryN).zfill(4) + ".png")
        if ( analysis ):
            outputFile.write("{0},{1},{2}\n".format(time/plotEveryN,dragCoeff,liftCoeff))
endTime = datetime.datetime.now()
deltaTime = endTime - startTime
if analysis:
    timeFile = open("{0}_time".format(prefix),"w")
    timeFile.write("{0}".format(deltaTime.total_seconds()));
    timeFile.write("\n");
    timeFile.close()
    outputFile.close()
os.chdir(workingFolder)
print 'End of calculation with {0} collision with Re={1} and size {2}.\n Elapsed time: {3}'.format(collStr,Re,size,deltaTime.total_seconds())
