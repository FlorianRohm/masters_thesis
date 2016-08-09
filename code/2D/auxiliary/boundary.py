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

from numpy import zeros

def cumulantBoundary(rho, u):
    ux = u[0, :]
    uy = u[1, :]
    fin = zeros((9,ux.size))

    uyq = uy**2
    uxq = ux**2
    zux = 2*ux
    zuy = 2*uy

    uxqqmux = (uxq - ux)
    uxqqpux = (uxq + ux)
    uyqqmuy = (uyq - uy)
    uyqqpuy = (uyq + uy)

    emuxqq = (1 - uxq)
    emuyqq = (1 - uyq)

    zuyme = (zuy - 1)
    zuype = (zuy + 1)


    K_20_p = 1./3.*rho
    K_02_p = K_20_p
    K_11_p = 0
    #back transformation
    c_22 = 2*K_11_p*K_11_p/rho + K_20_p*K_02_p/rho

    c_0_0 = (rho*uxqqmux  + K_20_p) * 0.5
    c_1_0 =  rho*emuxqq   - K_20_p
    c_2_0 = (rho*uxqqpux  + K_20_p) * 0.5
    c_0_1 = ( K_11_p*(zux - 1) ) * 0.5
    c_1_1 =  - K_11_p*zux
    c_2_1 = ( K_11_p*(zux + 1) ) * 0.5
    c_0_2 = (K_02_p*uxqqmux + c_22) * 0.5
    c_1_2 =  K_02_p*emuxqq  - c_22
    c_2_2 = (K_02_p*uxqqpux + c_22) * 0.5

    fin[7] = (c_0_0*uyqqmuy + c_0_1*zuyme + c_0_2) * 0.5
    fin[3] =  c_0_0*emuyqq  - c_0_1*zuy   - c_0_2
    fin[6] = (c_0_0*uyqqpuy + c_0_1*zuype + c_0_2) * 0.5
    fin[4] = (c_1_0*uyqqmuy + c_1_1*zuyme + c_1_2) * 0.5
    fin[0] =  c_1_0*emuyqq  - c_1_1*zuy   - c_1_2
    fin[2] = (c_1_0*uyqqpuy + c_1_1*zuype + c_1_2) * 0.5
    fin[8] = (c_2_0*uyqqmuy + c_2_1*zuyme + c_2_2) * 0.5
    fin[1] =  c_2_0*emuyqq  - c_2_1*zuy   - c_2_2
    fin[5] = (c_2_0*uyqqpuy + c_2_1*zuype + c_2_2) * 0.5

    return fin
