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

def distributionsFromCentralMoments_min (u, c_00, c_11, c_20, c_02, c_22):
    ux = u[0, :, :]
    uy = u[1, :, :]
    (nx,ny) = ux.shape
    f = zeros((3,3,nx,ny))

    uxq = ux**2
    zux = 2*ux
    uxqmux = (uxq - ux)
    uxqpux = (uxq + ux)
    emuxq = (1 - uxq)

    c_0_0 = (c_00*uxqmux  + c_20) * 0.5
    c_1_0 =  c_00*emuxq   - c_20
    c_2_0 = (c_00*uxqpux  + c_20) * 0.5
    c_0_1 = ( c_11*(zux - 1) ) * 0.5
    c_1_1 =  - c_11*zux
    c_2_1 = ( c_11*(zux + 1) ) * 0.5
    c_0_2 = (c_02*uxqmux + c_22) * 0.5
    c_1_2 =  c_02*emuxq  - c_22
    c_2_2 = (c_02*uxqpux + c_22) * 0.5

    uyq = uy**2
    zuy = 2*uy

    uyqmuy = (uyq - uy)
    zuyme = (zuy - 1)
    emuyq = (1 - uyq)
    zuype = (zuy + 1)
    uyqpuy = (uyq + uy)

    f[0, 0] = (c_0_0*uyqmuy + c_0_1*zuyme + c_0_2) * 0.5
    f[0, 1] =  c_0_0*emuyq  - c_0_1*zuy   - c_0_2
    f[0, 2] = (c_0_0*uyqpuy + c_0_1*zuype + c_0_2) * 0.5
    f[1, 0] = (c_1_0*uyqmuy + c_1_1*zuyme + c_1_2) * 0.5
    f[1, 1] =  c_1_0*emuyq  - c_1_1*zuy   - c_1_2
    f[1, 2] = (c_1_0*uyqpuy + c_1_1*zuype + c_1_2) * 0.5
    f[2, 0] = (c_2_0*uyqmuy + c_2_1*zuyme + c_2_2) * 0.5
    f[2, 1] =  c_2_0*emuyq  - c_2_1*zuy   - c_2_2
    f[2, 2] = (c_2_0*uyqpuy + c_2_1*zuype + c_2_2) * 0.5

    return array((f[1,1], f[2,1], f[1,2], f[0,1], f[1,0], f[2,2], f[0,2], f[0,0], f[2,0]))
