from cProfile import label
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
plt.style.use('seaborn-pastel')

fig = plt.figure()
ax = plt.axes()
line, = ax.plot([],[],'.')
f = open('data/energy.d', 'r')
f1 = open('data/gEnergy.d', 'r')
f2 = open('data/kEnergy.d', 'r')
dat=np.fromfile(f,dtype='float64', sep=" ").reshape(-1,2)
dat1=np.fromfile(f1,dtype='float64', sep=" ").reshape(-1,2)
dat2=np.fromfile(f2,dtype='float64', sep=" ").reshape(-1,2)

ax.plot(dat[:,0],dat[:,1], label = "$E$")
ax.plot(dat1[:,0],dat1[:,1], label ="$E_g$")
ax.plot(dat2[:,0],dat2[:,1], label = "$E_k$")
plt.legend()
plt.savefig('plots/energyPlot.png')
plt.xlabel("$time (s)$")
plt.ylabel("$Energy (j)$")
plt.show()

