from ctypes import *
from mpl_toolkits import mplot3d
from matplotlib import pyplot
tools = cdll.LoadLibrary('./tools.so')

ArrayType = c_double * 1600000

x = ArrayType()
y = ArrayType()
z = ArrayType()

ntot  = c_int(0)

filename = "/data2/mukherjeed2/starsmasher_nbody6_collisions/gam_updated/N_8000_S_0.8/coll_0/out0000.sph"

fn = filename.ljust(255)

tools.getxyzfromfile_(fn, byref(ntot), x, y, z)
print ntot

x_py = []
y_py = []
z_py = []

for i in xrange(0,(ntot.value-1)):
    x_py.append(x[i])
    y_py.append(y[i])
    z_py.append(z[i])

fig = pyplot.figure()
ax = fig.add_subplot(111, projection='3d')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

ax.scatter3D(x_py,y_py,z_py,s=0.5)
for angle in range(0, 360):
    ax.view_init(30, angle)
    pyplot.draw()
    pyplot.pause(.00001)
#pyplot.show()
