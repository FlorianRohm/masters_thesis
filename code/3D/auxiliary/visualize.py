from VTKWrapper import saveToVTK

def visualize(rho, u, time, saveEveryN, skipFirstN, grid, prefix):
    if ( (time % saveEveryN == 0) & (time > skipFirstN) ):

        saveNumber = str(time/saveEveryN).zfill(4)

        saveToVTK(rho, u, prefix, saveNumber, grid)
