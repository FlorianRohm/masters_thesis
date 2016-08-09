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
from matplotlib import cm, pyplot
from VTKWrapper import saveToVTK

def visualize(u, rho, timestep, plotEveryN, skipFirstN, liveUpdate, saveVTK, savePlot, ax, fig, grid, prefix):
    # Visualization
    if ( (timestep % plotEveryN == 0) & (liveUpdate | saveVTK | savePlot) & (timestep > skipFirstN) ):

        if ( liveUpdate | savePlot ):
            ax.clear()
            ax.imshow(sqrt(u[0]**2+u[1]**2).transpose(),  cmap=cm.afmhot, vmin=0., vmax=0.15)
            ax.set_title('velocity norm')

        if ( liveUpdate ):
            pyplot.draw()
        if ( saveVTK ):
            # convert velocity and density to 3d arrays
            printVel = reshape(u, (2, nx, ny, 1))
            printRho = reshape(rho, (nx, ny, 1))

            velocity = (printVel[0, :, :, :], printVel[1, :, :, :], velocityZ)
            saveNumber = str(timestep/plotEveryN).zfill(4)

            saveToVTK(velocity, printRho, prefix, saveNumber, grid)
        if ( savePlot ):
            pyplot.savefig(prefix + "." + str(timestep/plotEveryN).zfill(4) + ".png")
