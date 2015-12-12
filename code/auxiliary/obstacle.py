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
