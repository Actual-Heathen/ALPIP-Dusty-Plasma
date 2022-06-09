import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

N = 100
fig=plt.figure()
ax = plt.axes(xlim = (0,N), ylim = (0,N))

def animate(i):
    f = open('../data/density%05d.d' %(i), 'r')
    dat=np.fromfile(f, dtype='float64', sep= ' ').reshape(N,N)
    im = ax.imshow(dat)
    #plt.colormesh(dat)
    #plt.colorbar()
    f.close()

anim =  animation.FuncAnimation(fig,animate, frames = 200, interval = 20)

anim.save('density.gif',writer = 'imagemagick')