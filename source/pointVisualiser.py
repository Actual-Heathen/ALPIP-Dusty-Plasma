import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
plt.style.use('seaborn-pastel')

fig = plt.figure()
ax = plt.axes(xlim = (0,1/(100000)),ylim = (0,1/(100000)))
line, = ax.plot([],[],'.')

def init():
    line.set_data([],[])
    return line,

def animate(i):
    f = open('../data/points%05d.d' %(i), 'r')
    #for i in range (10):
    dat=np.fromfile(f,dtype='float64', sep=" ").reshape(-1,2)
    #print(dat)
        #x = np.linspace(0,4,1000)
        #y = np.sin(2*np.pi * (x-0.01*i)
    line.set_data(dat[:,0],dat[:,1])
    return line,

anim = FuncAnimation(fig,animate, init_func=init, frames = 200, interval = 20, blit = True)

anim.save('points.gif',writer = 'imagemagick')
