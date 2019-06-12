from pyfits import getdata
import matplotlib.pyplot as plt

mapname = '55.fits'

a = getdata(mapname)
plt.imshow(a,cmap = 'gray',vmin = 0,vmax = 30)
plt.axis('off')
plt.savefig('map.png')
