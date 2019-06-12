import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(4,18,1000)
y = np.ones(1000)
y = y*1.375
plt.plot(x,np.sin(x)+0.045*x,linestyle = ':')
plt.plot(x,np.sin(x)+0.045*x+0.3*np.sin(8.45983*x))
plt.plot(x,y,linestyle = '--',color = 'black')
plt.ylabel('Density $\\rho$')
plt.xlabel('$x$')
plt.xticks([])
plt.yticks([])
plt.savefig('bottomup.png')
