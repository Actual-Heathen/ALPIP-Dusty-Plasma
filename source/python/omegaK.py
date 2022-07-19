import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

k = np.linspace(-2,2, 1000)
w = np.zeros((2,1000),dtype=complex)
for ik in range (1000):
    w[:,ik]=np.roots([1,0,-16*k[ik]**2+1])

plt.plot(k,w[0,:].real)
plt.plot(k,w[0,:].imag)

plt.show()
plt.savefig('plots/omegaK.png')