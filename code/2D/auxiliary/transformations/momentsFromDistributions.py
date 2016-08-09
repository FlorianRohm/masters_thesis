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

def momentsFromDistributions (u, fin):
    ux = u[0, :, :]
    (nx,ny) = ux.shape
    f = zeros((3,3,nx,ny))

    f[1,1] = fin[0]
    f[2,1] = fin[1]
    f[1,2] = fin[2]
    f[0,1] = fin[3]
    f[1,0] = fin[4]
    f[2,2] = fin[5]
    f[0,2] = fin[6]
    f[0,0] = fin[7]
    f[2,0] = fin[8]
    c_0_0 =   f[0,0]    + f[0,1]     + f[0,2]
    c_0_1 = - f[0,0]    - f[0,1]     + f[0,2]
    c_0_2 =   f[0,0]    + f[0,1]     + f[0,2]
    c_1_0 =   f[1,0]    + f[1,1]     + f[1,2]
    c_1_1 = - f[1,0]    - f[1,1]     + f[1,2]
    c_1_2 =   f[1,0]    + f[1,1]     + f[1,2]
    c_2_0 =   f[2,0]    + f[2,1]     + f[2,2]
    c_2_1 = - f[2,0]    - f[2,1]     + f[2,2]
    c_2_2 =   f[2,0]    + f[2,1]     + f[2,2]
    c_00 =   c_0_0  + c_1_0  + c_2_0
    c_10 = - c_0_0  - c_1_0  + c_2_0
    c_20 =   c_0_0  + c_1_0  + c_2_0
    c_01 =   c_0_1  + c_1_1  + c_2_1
    c_11 = - c_0_1  - c_1_1  + c_2_1
    c_21 =   c_0_1  + c_1_1  + c_2_1
    c_02 =   c_0_2  + c_1_2  + c_2_2
    c_12 = - c_0_2  - c_1_2  + c_2_2
    c_22 =   c_0_2  + c_1_2  + c_2_2

    return (c_00, c_10, c_01, c_11, c_20, c_02, c_21, c_12, c_22)
