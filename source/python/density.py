import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

N = 4
N2 = N*N
fr = (20000)-2
fig=plt.figure()
ax = plt.axes(xlim = (0,N-1), ylim = (0,N-1))
d = open('data/densityF.d', 'r')
dat=np.fromfile(d, dtype='float64', sep= ' ', count=N2).reshape(N,N)
im = ax.imshow(dat)

def animate(i):
    print(i+1)
    dat=np.fromfile(d, dtype='float64', sep= ' ', count=N2).reshape(N,N)
    im.set_data(dat)
    #plt.colormesh(dat)

anim =  animation.FuncAnimation(fig,animate, frames = fr, interval =20)

anim.save('plots/density.gif',writer = 'imagemagick')
