#!/usr/bin/python
# 3D Lattice Boltzmann Code
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

for i in range(3):
    for j in range(3):
        print \
        "\
    k_{0}{1}_0 = f[{0},{1},0]              + f[{0},{1},1]           + f[{0},{1},2]          \n\
    k_{0}{1}_1 = f[{0},{1},0]*(-1 - uz)    + f[{0},{1},1]*(- uz)    + f[{0},{1},2]*(1 - uz) \n\
    k_{0}{1}_2 = f[{0},{1},0]*(-1 - uz)**2 + f[{0},{1},1]*(- uz)**2 + f[{0},{1},2]*(1 - uz)**2".format(i,j)

for i in range(3):
    for c in range(3):
        print \
        "\
    k_{0}_0{1} = k_{0}0_{1}              + k_{0}1_{1}           + k_{0}2_{1}          \n\
    k_{0}_1{1} = k_{0}0_{1}*(-1 - uy)    + k_{0}1_{1}*(- uy)    + k_{0}2_{1}*(1 - uy) \n\
    k_{0}_2{1} = k_{0}0_{1}*(-1 - uy)**2 + k_{0}1_{1}*(- uy)**2 + k_{0}2_{1}*(1 - uy)**2".format(i,c)

for b in range(3):
    for c in range(3):
        print \
        "\
    k_0{0}{1} = k_0_{0}{1}              + k_1_{0}{1}           + k_2_{0}{1}          \n\
    k_1{0}{1} = k_0_{0}{1}*(-1 - ux)    + k_1_{0}{1}*(- ux)    + k_2_{0}{1}*(1 - ux) \n\
    k_2{0}{1} = k_0_{0}{1}*(-1 - ux)**2 + k_1_{0}{1}*(- ux)**2 + k_2_{0}{1}*(1 - ux)**2".format(b,c)
