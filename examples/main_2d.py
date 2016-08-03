import sys, os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "../"))
from src.bobby import Bobby


def label(i, j, ni, nj):
    return i*(nj) + j

def write_2d_restart(ndim, nsize, nelem, x, y):
    f = open("restart.in", "w")
    f.write("%d\n"%ndim)
    f.write("%d\n"%nsize)
    f.write("%d\n"%nelem)
    x = x.ravel()
    y = y.ravel()
    for i in range(nsize):
        f.write("%d %.14e %.14e\n"%(i, x[i], y[i]))

    ne = int(np.sqrt(nelem))
    ns = int(np.sqrt(nsize))
    for i in range(ne):
        for j in range(ne):
            f.write("%d %d %d %d %d\n"%(label(i, j, ne, ne), label(i, j, ns, ns), label(i+1, j, ns, ns) ,label(i+1, j+1, ns, ns),label(i, j+1, ns, ns)))
            
    for i in range(nsize):
        if i == 0:
            f.write("%d %d %d\n"%(i, 0, nsize-1))
        elif i == nsize-1:
            f.write("%d %d %d\n"%(i, 0, 0))
        else:
            f.write("%d %d %d\n"%(i, -100, 0))
    f.close()

ndim = 2
nsize = 50
nelem = nsize - 1
xx = np.linspace(-1, 1, nsize)
yy = np.linspace(-1, 1, nsize)

x, y = np.meshgrid(xx, yy)
nsize = nsize**2
nelem = nelem**2
write_2d_restart(ndim, nsize, nelem, x, y)
