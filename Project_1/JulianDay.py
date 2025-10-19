import numpy as np

def JD(D,M,Y):
    a = np.floor((14 - M)/12)
    y = Y + 4800 - a
    m = M + 12*a -3

    JD = D + np.floor((153*m +2)/5 + 365*y +np.floor(y/4) - np.floor(y/100)) + np.floor(y/400) - 32045

    return JD