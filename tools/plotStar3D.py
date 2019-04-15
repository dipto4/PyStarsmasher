from ctypes import *
from mpl_toolkits import mplot3d
from matplotlib import pyplot
from optparse import OptionParser
import numpy as np

parser = OptionParser()
parser.add_option('-f','--filename',action='store',dest='filename')

(option,args) = parser.parse_args()

filename = str(option.filename)

tools = cdll.LoadLibrary('./tools.so')

ArrayType = c_double * 1600000

x = ArrayType()
y = ArrayType()
z = ArrayType()

ntot  = c_int(0)

#filename = "/data2/mukherjeed2/starsmasher_nbody6_collisions/gam_updated/N_8000_S_0.8/coll_0/out0000.sph"

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
#ax.axis("equal")
ax.set_aspect('equal')
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax.set_axis_off()
#test
max_range = np.array([max(x_py)-min(x_py), max(y_py)-min(y_py), max(z_py)-min(z_py)]).max() / 2.0

mid_x = (max(x_py)+min(x_py)) * 0.5
mid_y = (max(y_py)+min(y_py)) * 0.5
mid_z = (max(z_py)+min(z_py)) * 0.5
ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)

ax.scatter3D(x_py,y_py,z_py,s=0.5)
ax.view_init(30,30)
#for angle in range(0, 360):
#    ax.view_init(30, angle)
#    pyplot.draw()
#    pyplot.pause(.00001)
pyplot.show()
