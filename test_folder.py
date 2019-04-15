from PyStarsmasher import Starsmasher

folder = 'coll_0'

code = Starsmasher()

x = code.getLastOutFile(folder)

full_path = folder + '/' + x
boundMass = code.getBoundMass(full_path)
print(boundMass)
