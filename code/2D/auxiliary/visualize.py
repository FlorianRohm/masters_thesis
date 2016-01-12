from numpy import *
from matplotlib import cm, pyplot


def visualize(u, rho, timestep, plotEveryN, skipFirstN, liveUpdate, saveVTK, savePlot, ax, fig, prefix):
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
