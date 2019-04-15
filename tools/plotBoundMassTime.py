from ctypes import *
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-f','--folder',action='store',dest='folder')

(option, args) = parser.parse_args()

folder = str(option.folder)


tools = cdll.LoadLibrary('./tools.so')


Fboundmass = c_double(0.0)
for i in xrange(0,4001,10):
    filename = folder+'/out'+str(i).zfill(4)+'.sph'


    Pfilename = filename.ljust(255)
    tools.calculateboundmass_(Pfilename,byref(Fboundmass))
    print(Fboundmass.value)
    print i



