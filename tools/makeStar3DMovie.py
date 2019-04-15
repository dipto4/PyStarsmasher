from ctypes import *
from mpl_toolkits import mplot3d
from matplotlib import pyplot
from optparse import OptionParser
import numpy as np

def get_options():
    parser = OptionParser()
    parser.add_option('-f','--folder',action='store',dest='folder')
    parser.add_option('-s','--start',action='store',dest='start',default=0)
    parser.add_option('-e','--end',action='store',dest='end',default=4000)
    parser.add_option('-i','--interval',action='store',dest='interval',default=1)
    parser.add_option('-n','--numberofparticles',action='store',dest='num')

    (option,args) = parser.parse_args()
    folder = str(option.folder)
    start = int(option.start)
    end = int(option.end)
    interval = int(option.interaval)
    num = int(option.num)
    return folder,start,end,interval,num

def get_data_and_plot(folder,start,end,interval,num):

    count = 0
    tools = cdll.LoadLibrary('./tools.so')

    for i in xrange(start,(end+1),interval):

        ArrayType = c_double * 1600000

        x = ArrayType()
        y = ArrayType()
        z = ArrayType()

        ntot  = c_int(0)

        filename = folder+'/out'+str(i).zfill(4)+'.sph'
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

        if(count == 0):
            max_range = np.array([max(x_py)-min(x_py), max(y_py)-min(y_py), max(z_py)-min(z_py)]).max() / 2.0

            mid_x = (max(x_py)+min(x_py)) * 0.5
            mid_y = (max(y_py)+min(y_py)) * 0.5
            mid_z = (max(z_py)+min(z_py)) * 0.5

        plot(x_py,y_py,z_py,max_range,mid_x,mid_y,mid_z,num,count)
        count+=1

def plot(x_py,y_py,z_py,max_range,mid_x,mid_y,mid_z,num,count):

    fig = pyplot.figure()

    ax = fig.add_subplot(111, projection='3d')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
#ax.axis("equal")
    ax.set_aspect('equal')

#test
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    ax.scatter3D(x_py[0:(num-1)],y_py[0:(num-1)],z_py[0:(num-1)],s=0.5,c='blue')
    ax.scatter3D(x_py[num:],y_py[num:],z_py[num:],s=0.5,c='green')
    ax.view_init(30,30)


    filename = 'star3d_'+str(count).zfill(4)+'.png'
    pyplot.savefig(filename,dpi=600)
    pyplot.clf()


if __name__ == '__main__':
    folder,start,end,interval,num = get_options()
    get_data_and_plot(folder,start,end,interval,num)
