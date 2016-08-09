# 2D Lattice Boltzmann Code
# Copyright (C) 2016  Florian Rohm
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from numpy import *
from LBMHelpers import sumPopulations

# returns the boundaries of an obstacle as seen from particle distributions
def obstacleAttack(obstacle):
    (nx, ny) = obstacle.shape
    obstacleBound = zeros((9, nx, ny), dtype=bool)

    obstacleBound[1, :, :] = logical_xor(obstacle, logical_and(obstacle, roll(obstacle,  1, axis=0)))
    obstacleBound[2, :, :] = logical_xor(obstacle, logical_and(obstacle, roll(obstacle,  1, axis=1)))
    obstacleBound[3, :, :] = logical_xor(obstacle, logical_and(obstacle, roll(obstacle, -1, axis=0)))
    obstacleBound[4, :, :] = logical_xor(obstacle, logical_and(obstacle, roll(obstacle, -1, axis=1)))

    obstacleBound[5, :, :] = logical_xor(obstacle, logical_and(obstacle, roll(roll(obstacle,  1, axis=0),  1, axis=1)))
    obstacleBound[6, :, :] = logical_xor(obstacle, logical_and(obstacle, roll(roll(obstacle, -1, axis=0),  1, axis=1)))
    obstacleBound[7, :, :] = logical_xor(obstacle, logical_and(obstacle, roll(roll(obstacle, -1, axis=0), -1, axis=1)))
    obstacleBound[8, :, :] = logical_xor(obstacle, logical_and(obstacle, roll(roll(obstacle,  1, axis=0), -1, axis=1)))

    dragBounds = logical_or(obstacleBound[1, :, :],
                            logical_or(obstacleBound[5, :, :],
                                       logical_or(obstacleBound[8, :, :],
                                                  logical_or(obstacleBound[3, :, :],
                                                             logical_or(obstacleBound[6, :, :], obstacleBound[7, :, :])))))
    liftBounds = logical_or(obstacleBound[2, :, :],
                            logical_or(obstacleBound[6, :, :],
                                       logical_or(obstacleBound[5, :, :],
                                                  logical_or(obstacleBound[4, :, :],
                                                             logical_or(obstacleBound[7, :, :], obstacleBound[8, :, :])))))
    completeBound = logical_or(dragBounds, liftBounds)
    return (obstacleBound, dragBounds, liftBounds, completeBound)


def drag(scaledFin, bound):
    pos = sumPopulations(scaledFin[1, bound[1, :, :]]) + \
        sumPopulations(scaledFin[5, bound[5, :, :]]) + \
        sumPopulations(scaledFin[8, bound[8, :, :]])

    neg = sumPopulations(scaledFin[3, bound[3, :, :]]) + \
        sumPopulations(scaledFin[6, bound[6, :, :]]) + \
        sumPopulations(scaledFin[7, bound[7, :, :]])

    return pos-neg


def lift(scaledFin, bound):
    pos = sumPopulations(scaledFin[2, bound[2, :, :]]) + \
        sumPopulations(scaledFin[6, bound[6, :, :]]) + \
        sumPopulations(scaledFin[5, bound[5, :, :]])
    neg = sumPopulations(scaledFin[4, bound[4, :, :]]) + \
        sumPopulations(scaledFin[7, bound[7, :, :]]) + \
        sumPopulations(scaledFin[8, bound[8, :, :]])
    return pos-neg
