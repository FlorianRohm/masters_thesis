#!/usr/bin/python
# 6   2   5
#   \ | /
# 3 - 0 - 1
#   / | \
# 7   4   8

# 0,2   1,2   2,2
#   \    |   /
# 0,1   1,1   2,1
#   /    |   \
# 0,0   1,0   2,0

import os
import sys

filename = os.path.dirname(os.path.realpath(__file__)) + "/../transformations/distributionsFromCentralMoments.py"
fileObject = open(filename , "w")

orig_stdout = sys.stdout
sys.stdout = fileObject

print "#Automatically generated, do not change\n"
print "from numpy import *"
print ""
print "def distributionsFromCentralMoments (u, c_00, c_10, c_01, c_11, c_20, c_02, c_21, c_12, c_22):"
print "    ux = u[0, :, :]"
print "    uy = u[1, :, :]"
print "    (nx,ny) = ux.shape"
print "    f = zeros((3,3,nx,ny))"

for b in range(3):
    print \
    "\
    c_0_{0} = (c_0{0}*(ux**2 - ux) + c_1{0}*(2*ux - 1) + c_2{0}) * 0.5 \n\
    c_1_{0} =  c_0{0}*(1 - ux**2)  - c_1{0}*2*ux       - c_2{0}        \n\
    c_2_{0} = (c_0{0}*(ux**2 + ux) + c_1{0}*(2*ux + 1) + c_2{0}) * 0.5".format(b)

for i in range(3):
    print \
    "\
    f[{0}, 0] = (c_{0}_0*(uy**2 - uy) + c_{0}_1*(2*uy - 1) + c_{0}_2) * 0.5 \n\
    f[{0}, 1] =  c_{0}_0*(1 - uy**2)  - c_{0}_1*2*uy       - c_{0}_2        \n\
    f[{0}, 2] = (c_{0}_0*(uy**2 + uy) + c_{0}_1*(2*uy + 1) + c_{0}_2) * 0.5".format(i)

print ""
print "    return array((f[1,1], f[2,1], f[1,2], f[0,1], f[1,0], f[2,2], f[0,2], f[0,0], f[2,0]))"

sys.stdout = orig_stdout
fileObject.close()
