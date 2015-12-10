from LBMHelpers import iLeft, iCentV, iRight, iTop, iCentH, iBot

def YuLeft(f, feq):
    f[iRight, :] = feq[iLeft, :] + (feq[iRight, :] - f[iLeft, :])
