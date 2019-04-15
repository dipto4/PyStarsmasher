from ctypes import byref, cdll, c_int

lib = cdll.LoadLibrary('./libstarsmasher.so')

lib.testprint_()
