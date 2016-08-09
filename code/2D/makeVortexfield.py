#!/usr/bin/python
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

from auxiliary.vortex import buildVortexField
import sys, getopt
import os

import numpy

try:
  opts, args = getopt.getopt(sys.argv[1:],"hs:n:",)
except getopt.GetoptError:
    print "parse error"
    print 'test.py -s <edge length> -n <number of vortex pairs in one row/column> '
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print 'test.py -s <edge length> -n <number of vortex pairs in one row/column> '
        sys.exit()
    elif opt in ("-s"):
        length = int(float(arg))
    elif opt in ("-n"):
        nrOfPairs = int(float(arg))

diameter = length/(4.*nrOfPairs)

vel = buildVortexField(length,length,nrOfPairs,diameter)

numpy.save("l{0}_nr{1}_dia{2}".format(length,nrOfPairs,diameter),vel)
