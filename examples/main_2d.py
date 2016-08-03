import sys, os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "../"))
from src.bobby import Bobby
ndim = 2
nsize = 30
nn = nsize
ne = nn - 1


def label(i, j, ni, nj):
    return j*(nj) + i

def faces():
    bottom_face = [[i, i+1] for i in range(ne)]
    top_face = [[nn*(nn-1) + i+1-1, nn*(nn-1) + i-1] for i in range(ne, 0, -1)]
    right_face = [[nn*i+nn-1, nn*(i+1)+nn-1] for i in range(ne)]
    left_face = [[nn*(i), nn*(i-1)] for i in range(ne, 0, -1)]
    faces = []
    for b in [bottom_face, right_face, top_face, left_face]:
        faces.extend(b)
    return faces
def write_2d_restart(ndim, nsize, nelem, nb, x, y):
    f = open("restart.in", "w")
    f.write("%d\n"%ndim)
    f.write("%d\n"%nsize)
    f.write("%d\n"%nelem)
    f.write("%d\n"%nb)
    x = x.ravel()
    y = y.ravel()
    for i in range(nsize):
        f.write("%d %.14e %.14e\n"%(i, x[i], y[i]))

    ne = int(np.sqrt(nelem))
    ns = int(np.sqrt(nsize))
    for i in range(ne):
        for j in range(ne):
            f.write("%d %d %d %d %d\n"%(label(i, j, ne, ne), label(i, j, ns, ns), label(i+1, j, ns, ns) ,label(i+1, j+1, ns, ns),label(i, j+1, ns, ns)))
          
    _faces = faces()
    for i in range(nb):
        f.write("%d %d %d 0\n"%(i, _faces[i][0], _faces[i][1]))
    f.close()

nelem = nsize - 1
xx = np.linspace(-.5, .5, nsize)
yy = np.linspace(-.5, .5, nsize)
nb = nelem*4
x, y = np.meshgrid(xx, yy)

xr = x.copy()
yr = y.copy()

theta = np.pi/4.0

x =  xr*np.cos(theta) - yr*np.sin(theta)
y =  xr*np.sin(theta) + yr*np.cos(theta) 

nsize = nsize**2
nelem = nelem**2
write_2d_restart(ndim, nsize, nelem, nb, x, y)
