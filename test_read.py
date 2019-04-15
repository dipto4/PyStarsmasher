from optparse import OptionParser
from PyStarsmasher import Starsmasher

test = Starsmasher()

parser = OptionParser()
parser.add_option('-f','--folder',action="store",dest="folder")

(option, args) = parser.parse_args()
folder = str(option.folder)


rco_name, collision_track = test.getDataFromNBODY6(folder)

print("-------RESULTS---------")
print(rco_name)
#print(collision_track)
print len(collision_track)
for c in collision_track:
    print c

#print(main_collnums)
#print(main_names)
