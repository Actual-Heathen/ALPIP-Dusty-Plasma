import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
plt.style.use('seaborn-pastel')

fig = plt.figure()
ax = plt.axes(xlim = (0,200),ylim = (0,10/(1)))
#line, = ax.plot([],[],'.')
f = open('../data/pointTime.d', 'r')
dat=np.fromfile(f,dtype='float64', sep=" ").reshape(-1,2)
ax.plot(dat[:,0],dat[:,1], '.')
plt.savefig('pointTime.png')
plt.show()

