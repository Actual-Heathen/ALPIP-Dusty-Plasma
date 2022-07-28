from cProfile import label
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
plt.style.use('seaborn-pastel')

fig = plt.figure()
ax = plt.axes()
line, = ax.plot([],[],'.')
f = open('data/rhoS.d', 'r')

dat=np.fromfile(f,dtype='float64', sep=" ").reshape(-1,2)

ax.plot(dat[:,0],dat[:,1], label = "$rho^2$")
plt.legend()
plt.savefig('plots/rhoS.png')
plt.xlabel("$time (s)$")
plt.ylabel("$Energy (j)$")
plt.show()

