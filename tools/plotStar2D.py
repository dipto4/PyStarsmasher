from ctypes import *
from mpl_toolkits import mplot3d
from matplotlib import pyplot
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-f','--file',action='store',dest='outfile')
(option, args) = parser.parse_args()

f = str(option.outfile)

tools = cdll.LoadLibrary('./tools.so')

ArrayType = c_double * 1600000

x = ArrayType()
y = ArrayType()
z = ArrayType()

ntot  = c_int(0)


fn = f.ljust(255)
print(fn)
tools.getxyzfromfile_(fn, byref(ntot), x, y, z)
print ntot

x_py = []
y_py = []
z_py = []

for i in xrange(0,(ntot.value-1)):
    x_py.append(x[i])
    y_py.append(y[i])
    z_py.append(z[i])


pyplot.xlabel('x [RSun]',fontweight='bold')
pyplot.ylabel('y [RSun]',fontweight='bold')

pyplot.scatter(x_py,y_py,s=0.5)

pyplot.show()
#for angle in range(0, 360):
#    ax.view_init(30, angle)
#    pyplot.draw()
#    pyplot.pause(.00001)
#pyplot.show()


