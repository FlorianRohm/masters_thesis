#!/usr/bin/python

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
