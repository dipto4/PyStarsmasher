from ctypes import *
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-f','--file',action='store',dest='outfile')
(option, args) = parser.parse_args()

filename = str(option.outfile)



tools = cdll.LoadLibrary('./tools.so')


Pfilename = filename.ljust(255)
Fboundmass = c_double(0.0)
tools.calculateboundmass_(Pfilename,byref(Fboundmass))
print(Fboundmass.value)



