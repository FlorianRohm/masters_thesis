from pyevtk.hl import gridToVTK

def saveToVTK(rho, u, prefix, saveNumber, grid):
    name = "./" + prefix + "." + saveNumber
    reorderedU = (u[0], u[1], u[2])
    gridToVTK(name, grid[0], grid[1], grid[2],
              pointData = {'velocity': reorderedU,
                           'pressure': rho})
