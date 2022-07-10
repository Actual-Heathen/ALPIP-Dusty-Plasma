import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
plt.style.use('seaborn-pastel')

fig = plt.figure()
ax = plt.axes(xlim = (0,1/(1)),ylim = (0,1/(1)))
line, = ax.plot([],[],'.')
f = open('data/pointsF.d', 'r')
particles = 2500

def init():
    line.set_data([],[])
    return line,

def animate(i):

    #for i in range (10):
    dat=np.fromfile(f,dtype='float64', sep=" ", count = (particles*2)).reshape(-1,2)
    #print(dat)
        #x = np.linspace(0,4,1000)
        #y = np.sin(2*np.pi * (x-0.01*i)
    line.set_data(dat[:,0],dat[:,1])
    return line,

anim = FuncAnimation(fig,animate, init_func=init, frames = 500, interval = 20, blit = True)

anim.save('plots/points.gif',writer = 'imagemagick')
