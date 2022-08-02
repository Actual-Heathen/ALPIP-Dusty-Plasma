import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as colors
import matplotlib as mpl

d = open('data/density.d', 'r')
f = open('data/fTransform.d', 'r')
a = open('data/aFTransform.d', 'r')
p = open('data/psi.d', 'r')

N = 100
N2 = N*N
fr = (1000)-2

fig = plt.figure()

def update_anim(it):
    
    fig.clf() #clear the figure
    
    ax=fig.subplots(2, 2) #add subplots

    fig.tight_layout() #reduce spacing around the figure

    #for i in range(2):
    #    for j in range(2):
    #        ax[i,j].cla()  

                    
    Z1=np.fromfile(d, dtype='float64', sep= ' ', count=N2).reshape(N,N)
    Z2=np.fromfile(f, dtype='float64', sep= ' ', count=N2).reshape(N,N)
    Z3=np.fromfile(a, dtype='float64', sep= ' ', count=N2).reshape(N,N)
    Z4=np.fromfile(p, dtype='float64', sep= ' ', count=N2).reshape(N,N)

    im1=ax[0,0].imshow(Z1 ,aspect='equal', origin='lower',cmap='plasma', norm=mpl.colors.SymLogNorm(linthresh=1, vmax=100,vmin=0))
    im2=ax[1,0].imshow(Z2 ,aspect='equal', origin='lower', cmap='Blues')
    im3=ax[0,1].imshow(Z3 ,aspect='equal', origin='lower', cmap='Reds')
    im4=ax[1,1].imshow(Z4 ,aspect='equal', origin='lower', cmap='Greens')
    #im4=ax[1,1].plot(x,np.sin(it*x))

    

    

    cb1=fig.colorbar(im1, ax=ax[0,0],extend='both')
    cb2=fig.colorbar(im2, ax=ax[1,0])
    cb3=fig.colorbar(im3, ax=ax[0,1])
    cb4=fig.colorbar(im4, ax=ax[1,1])
    print(it)


fig=plt.figure(figsize=(12,12))
anim=animation.FuncAnimation(fig,update_anim,frames=fr,interval=20)

anim.save('plots/multiDensity.gif')
plt.subplots_adjust(wspace=0.1,hspace=0.01, bottom=0,top=1)

plt.close()
anim