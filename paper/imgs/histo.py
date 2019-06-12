import numpy as np
import matplotlib.pyplot as plt
from pyfits import getdata
from scipy import signal, ndimage

rad = 2
a = getdata('10.fits')
s = np.size(a,0)
a2 = np.zeros((s,s))
a5 = a2
a10 = a2
ndimage.gaussian_filter(a,[2,2],output=a2)
ndimage.gaussian_filter(a,[5,5],output=a5)
ndimage.gaussian_filter(a,[10,10],output=a10)
a = a.ravel()
a2 = a2.ravel()
a5 = a5.ravel()
a10 = a10.ravel()

histo,edges = np.histogram(a,bins=500,range=(0,16), normed=True)
histo2 = np.histogram(a2,bins=500,range=(0,10), normed=True)[0]
histo5 = np.histogram(a5,bins=500,range=(0,10), normed=True)[0]
histo10 = np.histogram(a10,bins=500,range=(0,16), normed=True)[0]

plt.plot(edges[1:],histo, label='Unconvolved Map')
#plt.plot(edges[1:],histo2)
#plt.plot(edges[1:],histo5, label='Convolved Map')
plt.plot(edges[1:],histo10, label='Convolved Map')
plt.legend()
plt.xlabel('Magnification (linear)')
plt.ylabel('Probability')
plt.savefig('convhisto.png')
plt.show()
