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
#
from numpy import *

def centralMomentsFromDistributions_min (u, fin):
    ux = u[0, :, :]
    uy = u[1, :, :]
    (nx,ny) = ux.shape
    f = zeros((3,3,nx,ny))
    uq = uy**2
    emu = (1 - uy)
    emuq = emu**2
    memu = (-1 - uy)
    memuq = memu**2

    f[1,1] = fin[0]
    f[2,1] = fin[1]
    f[1,2] = fin[2]
    f[0,1] = fin[3]
    f[1,0] = fin[4]
    f[2,2] = fin[5]
    f[0,2] = fin[6]
    f[0,0] = fin[7]
    f[2,0] = fin[8]
    c_0_0 = f[0,0]       + f[0,1]     + f[0,2]
    c_0_1 = f[0,0]*memu  - f[0,1]*uy  + f[0,2]*emu
    c_0_2 = f[0,0]*memuq + f[0,1]*uq  + f[0,2]*emuq
    c_1_0 = f[1,0]       + f[1,1]     + f[1,2]
    c_1_1 = f[1,0]*memu  - f[1,1]*uy  + f[1,2]*emu
    c_1_2 = f[1,0]*memuq + f[1,1]*uq  + f[1,2]*emuq
    c_2_0 = f[2,0]       + f[2,1]     + f[2,2]
    c_2_1 = f[2,0]*memu  - f[2,1]*uy  + f[2,2]*emu
    c_2_2 = f[2,0]*memuq + f[2,1]*uq  + f[2,2]*emuq

    memux = (-1 - ux)
    emux = (1 - ux)

    c_20 = c_0_0*memux**2 + c_1_0*ux**2     + c_2_0*emux**2
    c_11 = c_0_1*memux    - c_1_1*ux        + c_2_1*emux
    c_02 = c_0_2              + c_1_2           + c_2_2

    return (c_11, c_20, c_02)
