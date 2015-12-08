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
# 5   2   8
#   \ | /
# 3 - 0 - 6
#   / | \
# 4   1   7
#

from numpy import *
from matplotlib import cm, pyplot
from VTKWrapper import saveToVTK
from ghiaResults import *
import os


###### Plot settings ############################################################

plotEveryN    = 10         # draw every plotEveryN'th cycle
skipFirstN    = 0       # do not process the first skipFirstN cycles
savePlot      = True      # save velocity norm and x velocity plot
liveUpdate    = True     # show the process of the simulation (slow)
saveVTK       = False       # save the vtk files
prefix        = 'ldc'      # naming prefix for saved files
outputFolder  = './out'    # folder to save the output to
workingFolder = os.getcwd()


###### Flow definition #########################################################
maxIterations = 1  # Total number of time iterations.
Re            = 100.0   # Reynolds number.re

# Number of Cells
nx = 200
ny = 200

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


###### Lattice Constants #######################################################

# lattice velocities
c = array([[0, 0], [1, 0], [0, 1], [-1, 0], [0, -1], [1, 1], [-1, 1], [-1, -1], [1, -1]])
# Lattice weights
t      = 1./36. * ones(q)
t[1:5] = 1./9.
t[0]   = 4./9.

# index array for noslip reflection
# inverts the d2q9 stencil for use in bounceback scenarios
noslip = [0, 3, 4, 1, 2, 7, 8, 5, 6]

# index arrays for different sides of the d2q9 stencil
iLeft   = arange(q)[asarray([ci[0] <  0 for ci in c])]
iCentV  = arange(q)[asarray([ci[0] == 0 for ci in c])]
iRight  = arange(q)[asarray([ci[0] >  0 for ci in c])]
iTop    = arange(q)[asarray([ci[1] >  0 for ci in c])]
iCentH  = arange(q)[asarray([ci[1] == 0 for ci in c])]
iBot    = arange(q)[asarray([ci[1] <  0 for ci in c])]


###### Function Definitions ####################################################

# Helper function for density computation.
def sumPopulations(fin):
    return sum(fin, axis = 0)

def collide(fin, omega, u):
    print "precollide"

    print amax(fin[0, :, :])

    print amax(fin[1, :, :])
    print amax(fin[2, :, :])
    print amax(fin[3, :, :])
    print amax(fin[4, :, :])

    print amax(fin[5, :, :])
    print amax(fin[6, :, :])
    print amax(fin[7, :, :])
    print amax(fin[8, :, :])

    print "\ncollide"
    # central moments
    ux = u[0, :, :]
    uy = u[1, :, :]

    print "\nvelocities"
    print amax(ux)
    print amax(uy)

    # c_{i beta}
    # c_{-1 beta}
    c__n1_0 = fin[7, :, :] + fin[4, :, :] + fin[8, :, :] # f-1-1 + f-10 + f-11
    c__n1_1 = (-1-uy)*fin[7, :, :] - uy*fin[4, :, :] + (1-uy)*fin[8, :, :]
    c__n1_2 = (-1-uy)*(-1-uy)*fin[7, :, :] + uy*uy*fin[4, :, :] + (1-uy)*(1-uy)*fin[8, :, :]

    # c_{0 beta}
    c__0_0 = fin[3, :, :] + fin[0, :, :] + fin[1, :, :] # f0-1 + f00 + f01
    c__0_1 = (-1-uy)*fin[3, :, :] - uy*fin[0, :, :] + (1-uy)*fin[1, :, :]
    c__0_2 = (-1-uy)*(-1-uy)*fin[3, :, :] + uy*uy*fin[0, :, :] + (1-uy)*(1-uy)*fin[1, :, :]

    # c_{1 beta}
    c__1_0 = fin[6, :, :] + fin[2, :, :] + fin[5, :, :] # f1-1 + f10 + f11
    c__1_1 = (-1-uy)*fin[6, :, :] - uy*fin[2, :, :] + (1-uy)*fin[5, :, :]
    c__1_2 = (-1-uy)*(-1-uy)*fin[6, :, :] + uy*uy*fin[2, :, :] + (1-uy)*(1-uy)*fin[5, :, :]

    # c{alpha beta}
    # c{alpha 0}
    c_00 = c__n1_0 + c__0_0 + c__1_0
    c_10 = (-1-ux)*c__n1_0 - ux*c__0_0 + (1-ux)*c__1_0
    c_20 = (-1-ux)*(-1-ux)*c__n1_0 + ux*ux*c__0_0 + (1-ux)*(1-ux)*c__1_0

    # c{alpha 1}
    c_01 = c__n1_1 + c__0_1 + c__1_1
    c_11 = (-1-ux)*c__n1_1 - ux*c__0_1 + (1-ux)*c__1_1
    c_21 = (-1-ux)*(-1-ux)*c__n1_1 + ux*ux*c__0_1 + (1-ux)*(1-ux)*c__1_1

    # c{alpha 2}
    c_02 = c__n1_2 + c__0_2 + c__1_2
    c_12 = (-1-ux)*c__n1_2 - ux*c__0_2 + (1-ux)*c__1_2
    c_22 = (-1-ux)*(-1-ux)*c__n1_2 + ux*ux*c__0_2 + (1-ux)*(1-ux)*c__1_2

    # only K22 differs from the central moments
    K_22 = c_22 - 2*c_11*c_11 - c_20*c_02

    print "\nnormalized cumulants"
    print amax(c_00)
    print amax(c_10)
    print amax(c_01)
    print amax(c_11)
    print amax(c_20)
    print amax(c_02)
    print amax(c_21)
    print amax(c_12)
    print amax(c_22)

    # collision
    c_00_post = c_00
    c_10_post = c_10
    c_01_post = c_01

    c_11_post = (1-omega)*c_11

    c_20_post = 1/3 + 0.5*(1-omega)*(c_20 - c_02)
    c_02_post = 1/3 - 0.5*(1-omega)*(c_20 - c_02)

    c_21_post = 0
    c_12_post = 0

    K_22_post = 0

    # Transformation to central moments
    c_22_post = K_22_post + 2*c_11_post*c_11_post + c_20_post*c_02_post

    # backward transformation
    c__n1_0_post = -0.5*(ux*(1 - ux))*c_00_post - 0.5*(1-2*ux)*c_10_post + 0.5*c_20_post
    c__n1_1_post = -0.5*(ux*(1 - ux))*c_01_post - 0.5*(1-2*ux)*c_11_post + 0.5*c_21_post
    c__n1_2_post = -0.5*(ux*(1 - ux))*c_02_post - 0.5*(1-2*ux)*c_12_post + 0.5*c_22_post

    c__0_0_post = (1-ux*ux)*c_00_post - 2*ux*c_10_post - c_20_post
    c__0_1_post = (1-ux*ux)*c_01_post - 2*ux*c_11_post - c_21_post
    c__0_2_post = (1-ux*ux)*c_02_post - 2*ux*c_12_post - c_22_post

    c__1_0_post = 0.5*(ux*(1 - ux))*c_00_post + 0.5*(1+2*ux)*c_10_post + 0.5*c_20_post
    c__1_1_post = 0.5*(ux*(1 - ux))*c_01_post + 0.5*(1+2*ux)*c_11_post + 0.5*c_21_post
    c__1_2_post = 0.5*(ux*(1 - ux))*c_02_post + 0.5*(1+2*ux)*c_12_post + 0.5*c_22_post

    # post collision distributions
    f__n1_n1_post = -0.5*(uy*(1 - uy))*c__n1_0_post - 0.5*(1-2*uy)*c__0_0_post + 0.5*c__1_0_post
    f__n1_0_post = -0.5*(uy*(1 - uy))*c__n1_1_post - 0.5*(1-2*uy)*c__0_1_post + 0.5*c__1_1_post
    f__n1_1_post = -0.5*(uy*(1 - uy))*c__n1_2_post - 0.5*(1-2*uy)*c__0_2_post + 0.5*c__1_2_post

    f__0_n1_post = (1-uy*uy)*c__n1_0_post - 2*uy*c__0_0_post - c__1_0_post
    f__0_0_post = (1-uy*uy)*c__n1_1_post - 2*uy*c__0_1_post - c__1_1_post
    f__0_1_post = (1-uy*uy)*c__n1_2_post - 2*uy*c__0_2_post - c__1_2_post

    f__1_n1_post = 0.5*(uy*(1 - uy))*c__n1_0_post + 0.5*(1+2*uy)*c__0_0_post + 0.5*c__1_0_post
    f__1_0_post = 0.5*(uy*(1 - uy))*c__n1_1_post + 0.5*(1+2*uy)*c__0_1_post + 0.5*c__1_1_post
    f__1_1_post = 0.5*(uy*(1 - uy))*c__n1_2_post + 0.5*(1+2*uy)*c__0_2_post + 0.5*c__1_2_post

    print "\npost collision cumulants"
    print amax(f__0_0_post)

    print amax(f__1_0_post)
    print amax(f__0_1_post)
    print amax(f__n1_0_post)
    print amax(f__0_n1_post)

    print amax(f__1_1_post)
    print amax(f__n1_1_post)
    print amax(f__n1_n1_post)
    print amax(f__1_n1_post)

    return array([f__0_0_post, f__1_0_post, f__0_1_post, f__n1_0_post, f__0_n1_post,   f__1_1_post, f__n1_1_post, f__n1_n1_post, f__1_n1_post])

# Equilibrium distribution function.
def equilibrium(rho, u):
    cu   = 3.0 * dot(c, u.transpose(1, 0, 2))
    usqr = 3./2.*(u[0]**2+u[1]**2)
    feq = zeros((q, nx, ny))
    for i in range(q):
        feq[i, :, :] = rho*t[i]*(1.+cu[i]+0.5*cu[i]**2-usqr)
        # feq[i, :, :] = t[i]*(rho+cu[i]+0.5*cu[i]**2-usqr)
    return feq


def clamp(val, minVal, maxVal):
    return maximum(minVal, minimum(val, maxVal))

###### Setup ##################################################################

# set up the walls
leftWall    = fromfunction(lambda   x, y: x == 0,  (nx, ny))
rightWall   = fromfunction(lambda   x, y: x == nxl,  (nx, ny))
bottomWall  = fromfunction(lambda   x, y: y == nyl,  (nx, ny))

wall    = logical_or(logical_or(leftWall, rightWall), bottomWall)

# initial velocity profile
#  < -
#   -
#  - >
maxVel                = zeros((2, nx, ny))
maxVel[0, :,  0 ]     =  uLB

# initial particle distributions
feq   = equilibrium(1.0, maxVel)
fin   = feq.copy()
fpost = feq.copy()

# interactive mode (execute code while showing figures)
if ( liveUpdate | savePlot ):
    pyplot.ion()
    f = pyplot.figure(figsize=(15, 7.5))
    subplot1 = pyplot.subplot2grid((2, 4), (0, 0), colspan=2, rowspan=2)
    subplot2 = pyplot.subplot2grid((2, 4), (0, 2), colspan=2, rowspan=2)

os.chdir(outputFolder)


###### Main time loop ##########################################################
for time in range(maxIterations):

    # Calculate macroscopic density ...
    rho = sumPopulations(fin)
    # ... and velocity
    u = dot(c.transpose(),  fin.transpose((1, 0, 2)))/rho
    # u = dot(c.transpose(),  fin.transpose((1, 0, 2)))

    feq = equilibrium(rho, u)

    # Collision step.
    fpost = collide(fin, omega, u)

    # Streaming step
    fin[0, :, :] = fpost[0, :, :]

    fin[1, 1:nxl,   :]     = fpost[1, 0:nxl-1,  :]
    fin[2,   :,   0:nyl-1] = fpost[2,   :,    1:nyl]
    fin[3, 0:nxl-1, :]     = fpost[3, 1:nxl,    :]
    fin[4,   :,   1:nyl]   = fpost[4,   :,    0:nyl-1]

    fin[5, 1:nxl,   0:nyl-1] = fpost[5, 0:nxl-1, 1:nyl]
    fin[6, 0:nxl-1, 0:nyl-1] = fpost[6, 1:nxl,   1:nyl]
    fin[7, 0:nxl-1, 1:nyl]   = fpost[7, 1:nxl,   0:nyl-1]
    fin[8, 1:nxl,   1:nyl]   = fpost[8, 0:nxl-1, 0:nyl-1]


    # impulse ramping
    factor = clamp(time, 100, 100)/100.
    vel = maxVel * factor
    # Left wall: compute density from known populations
    u[:, :, 0] = vel[:, :, 0]
    rho[:, 0] = sumPopulations(fin[iCentV, 0, :])+2.*sumPopulations(fin[iTop, 0, :])


    # complete the left wall treatement wrt Yu 2002
    fin[iBot, 0, :] = - feq[iTop, 0, :] + (feq[iBot, 0, :] + fin[iTop, 0, :])
    # fin[iBot, 0, :] = fin[iTop, 0, :] + 6 * dot(c, u.transpose(1, 0, 2)) c[iBot][0] * rho*t[iBot]

    # wall boundary handling
    for i in range(q):
        fin[i, wall] = fin[noslip[i], wall]




    # Visualization
    if ( (time % plotEveryN == 0) & (liveUpdate | saveVTK | savePlot) & (time > skipFirstN) ):
        if ( liveUpdate | savePlot ):
            subplot1.clear()
            subplot2.clear()
            subplot1.imshow(sqrt(u[0]**2+u[1]**2).transpose(),  cmap=cm.jet, vmin=0., vmax=0.05)
            subplot1.set_title('velocity norm')

            velNormX = u[0, nx/2, :]/uLB
            subplot2.plot(velNormX, axisYNorm, label="lbm result")
            subplot2.plot(xghia, yghia, 'go' , label="ghia et al")
            subplot2.set_title('x velocity in center column')

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
