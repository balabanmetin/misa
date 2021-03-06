import numpy as np


def jc69(p):
    if (p - np.finfo(float).eps < 0):
        return 0.0
    else:
        loc = 1 - (4 * p / 3)
        if (0 >= loc):
            return  5.0 # treat as missing data
        else:
            return -0.75*np.log(loc)