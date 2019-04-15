from matplotlib import pyplot
from matplotlib import ticker
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-f','--folder',action='store',dest='folder')

(option,args) = parser.parse_args()


folder = str(option.folder)


time = []
W = []
T = []
U = []
E= []

energyfile = folder+'/'+'energy0.sph'

fenergy = open(energyfile,'r')

for line in fenergy:
    params = line.split()
    time.append(float(params[0]))
    W.append(float(params[1]))
    T.append(float(params[2]))
    U.append(float(params[3]))
    E.append(float(params[4]))

fig, axs = pyplot.subplots(4, 1,sharex=True)
#axs[0].grid(True)
#pyplot.title('Change in Energy of System Over Time',fontweight='bold')

axs[0].set_title('Change in Energy of System Over Time', fontweight='bold')


axs[0].yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
axs[0].yaxis.set_major_locator(pyplot.MaxNLocator(5))

axs[0].plot(time,W)
axs[0].set_ylabel('W',fontweight='bold')
axs[1].set_ylabel('T',fontweight='bold')
axs[2].set_ylabel('U',fontweight='bold')
axs[3].set_ylabel('E',fontweight='bold')
axs[3].set_xlabel('Time [code units]', fontweight='bold')
axs[1].plot(time,T)
axs[2].plot(time,U)
axs[3].plot(time,E)
fig.subplots_adjust(hspace=0.1)
for ax in axs:
    #ax.grid(True)
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.yaxis.set_major_locator(pyplot.MaxNLocator(5))
    #ax.yaxis.set_label_coords(-0.2,1.02)

pyplot.show()

