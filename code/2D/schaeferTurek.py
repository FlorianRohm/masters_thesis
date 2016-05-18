#!/usr/bin/python
# Copyright (C) 2013 FlowKit Ltd, Lausanne, Switzerland
# E-mail contact: contact@flowkit.com
#
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License, either
# version 3 of the License, or (at your option) any later version.
#
# Extended Comments and bugfixed by Florian Rohm, 2015

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
from auxiliary.collide import BGKCollide, cumulantCollide_min
from auxiliary.stream import stream
from auxiliary.LBMHelpers import clamp, getMacroValues, sumPopulations, equilibrium, noslip, iLeft, iCentV, iRight, iTop, iCentH, iBot
from auxiliary.boundaryConditions import YuLeft
from auxiliary.obstacle import obstacleAttack, drag, lift
import sys, getopt
import os

###### Parsing ############################################################

Re = 100
size = 10
collisionFunction = BGKCollide
try:
  opts, args = getopt.getopt(sys.argv[1:],"hcr:s:",["re=","size="])
except getopt.GetoptError:
  print 'test.py -i <inputfile> -o <outputfile>'
  sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print 'test.py -r <reynolds number, default 100> -s <size of the sphere, divisible by 10, default 10>'
        sys.exit()
    elif opt in ("-r" "--re"):
        Re = int(float(arg))
    elif opt in ("-s", "--size"):
        size = int(float(arg))
    elif opt in ("-c", "--cumulant"):
        collisionFunction = cumulantCollide_min
        print "Using cumulant collision"
    else:
        print "Using srt collision"
print 'Reynolds number is', Re
print 'Size of sphere is', size
factor = size/10.


###### Plot settings ############################################################

plotEveryN    = 100         # draw every plotEveryN'th cycle
skipFirstN    = 0       # do not process the first skipFirstN cycles
savePlot      = False      # save velocity norm and x velocity plot
liveUpdate    = True      # show the process of the simulation (slow)
analysis      = True       # write out drag and lift
prefix        = 'schaeferTurek_cumulant_Re{0}_size{1}'.format(Re, size)      # naming prefix for saved files
outputFolder  = './out'    # folder to save the outputFile to
workingFolder = os.getcwd()


###### Flow definition #########################################################
maxIterations = 20000  # Total number of time iterations.

# Number of Cells
ny = int(round(41*factor)) + 2 # for boundary
nx = int(round(220*factor)) + 2 # for boundary

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
outputFile = open(prefix, 'w')

outputFile.write("timestep,drag,lift\n")

###### Setup ##################################################################

# cylindrical obstacle
obstacle       = fromfunction(lambda    x, y: (x-cx)**2+(y-cy)**2 < r**2,  (nx, ny))
obstacleBounds, dragBoundStencil, liftBoundStencil, completeBoundStencil = obstacleAttack(obstacle)

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

    # Right Wall: Produce zero pressure gradient for the outflow
    fin[iLeft, -1, :] = fin[iLeft, -2, :]

    # Calculate macroscopic density and velocity
    (rho, u) = getMacroValues(fin)

    # Left wall: compute density from known populations.
    u[:, 0, :] = vel[:, 0, :]
    rho[0, :] = 1./(1.-u[0, 0, :]) * (sumPopulations(fin[iCentV, 0, :])+2.*sumPopulations(fin[iLeft, 0, :]))

    feq[:,0:1,:] = equilibrium(rho[0:1,:], u[:,0:1,:])

    # complete the left wall treatement wrt Yu 2002
    fin[iRight, 0, :] = feq[iLeft, 0, :] + (feq[iRight, 0, :] - fin[iLeft, 0, :])

    # Collision step.
    fpost[:,1:nx-1,1:ny-1] = collisionFunction(fin[:,1:nx-1,1:ny-1], rho[1:nx-1,1:ny-1], u[:,1:nx-1,1:ny-1], omega )

    # Streaming step
    fin = stream(fpost)

    # Visualization
    if ( (time % plotEveryN == 0) & (analysis | liveUpdate  | savePlot) & (time > skipFirstN) ):
        # Here, distributions are streamed into the obstacle -> compute drag and lift
        scaling = 2/(sumPopulations(fin)*uLBAverage*uLBAverage*size)
        scaledFin = fin.copy()
        # just scale the populations which are used
        for i in range(9):
            scaledFin[i, completeBoundStencil] = scaledFin[i, completeBoundStencil]*scaling[completeBoundStencil]


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
            outputFile.write("{0},{1},{2}\n".format(time,dragCoeff,liftCoeff))

outputFile.close()
os.chdir(workingFolder)
