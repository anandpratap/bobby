import sys, os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "../"))
from src.bobby import Bobby

def write_1d_restart(ndim, nsize, nelem, x):
    f = open("restart.in", "w")
    f.write("%d\n"%ndim)
    f.write("%d\n"%nsize)
    f.write("%d\n"%nelem)
    for i in range(nsize):
        f.write("%d %.14e\n"%(i, x[i]))

    for el in range(nelem):
        f.write("%d %d %d\n"%(el, el, el+1))

    for i in range(nsize):
        if i == 0:
            f.write("%d %d %d\n"%(i, 0, nsize-1))
        elif i == nsize-1:
            f.write("%d %d %d\n"%(i, 0, 0))
        else:
            f.write("%d %d %d\n"%(i, -100, 0))
    f.close()

ndim = 1
nsize = 300
nelem = nsize - 1
x = np.linspace(0, 1, nsize)
write_1d_restart(ndim, nsize, nelem, x)

b = Bobby()
b.run()

u = b.get_solution(nsize*3)
x = b.get_mesh(nsize)

plt.figure()
plt.plot(x, u[::b.nvar])
plt.plot(x, u[1::b.nvar])
plt.show()
