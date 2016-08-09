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

from numpy import *

def centralMomentsFromCumulants( rhoInv , k_000, k_001, k_002, k_010, k_011, k_012, k_020, k_021, C_022,
         k_100, k_101, k_102, k_110, k_111, C_112, k_120, C_121, C_122,
         k_200, k_201, C_202, k_210, C_211, C_212, C_220, C_221, C_222):

    k_211 = C_211 + (k_200*k_011 + 2*k_110*k_101)*rhoInv
    k_121 = C_121 + (k_020*k_101 + 2*k_011*k_110)*rhoInv
    k_112 = C_112 + (k_002*k_110 + 2*k_101*k_011)*rhoInv

    k_220 = C_220 + (k_200*k_020 + 2*k_110**2)*rhoInv
    k_202 = C_202 + (k_002*k_200 + 2*k_101**2)*rhoInv
    k_022 = C_022 + (k_020*k_002 + 2*k_011**2)*rhoInv

    k_122 = C_122 + (k_002*k_120 + k_020*k_120 + 4*k_011*k_111 + 2*(k_101*k_021 + k_110*k_012))*rhoInv
    k_212 = C_212 + (k_002*k_210 + k_200*k_210 + 4*k_101*k_111 + 2*(k_011*k_201 + k_110*k_102))*rhoInv
    k_221 = C_221 + (k_020*k_201 + k_200*k_201 + 4*k_110*k_111 + 2*(k_011*k_210 + k_101*k_120))*rhoInv

    k_222 = C_222 + (4*k_111**2 + k_200*k_022 + k_020*k_202 + k_002*k_220  \
                     + 4*(k_011*k_211 + k_101*k_121 + k_110*k_112) \
                     + 2*(k_120*k_102 + k_210*k_012 + k_201*k_021))*rhoInv \
                  - (16*k_110*k_101*k_011 + 4*(k_020*k_101**2 + k_200*k_011**2 + k_002*k_110**2) \
                     + 2*k_200*k_020*k_002)*rhoInv**2

    return ( k_000, k_001, k_002, k_010, k_011, k_012, k_020, k_021, k_022,
             k_100, k_101, k_102, k_110, k_111, k_112, k_120, k_121, k_122,
             k_200, k_201, k_202, k_210, k_211, k_212, k_220, k_221, k_222)
