import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u

n = 1000
rscale = 10.0 * u.pc
rlim   = 300.0 * u.pc

xcenter = 500.0 * u.pc
ycenter = 400.0 * u.pc

R = np.random.rand(n) * (1.0 - np.exp(-rlim/rscale)) 
rvals = -1.0 * np.log(1.0 - R) * rscale 
phi = np.random.rand(n)*2.0*np.pi

x = np.cos(phi)*rvals # 
y = np.sin(phi)*rvals # in dwarf coordinates

# convert to 
x = x + xcenter
y = y + ycenter

plt.scatter(x.value,y.value,s=5,alpha=0.5)
plt.xlabel('x'); plt.ylabel('y')
plt.xlim(0,2000.0); plt.ylim(0.,1000.0)
plt.axis('equal')
circle = np.linspace(0.0,2.0*np.pi)
plt.plot(rlim.value * np.cos(circle) + xcenter.value, rlim.value * np.sin(circle)+ycenter.value, color='red')
plt.savefig('xy_scatter.png')
plt.close()

nbins = 25
hist_vals = plt.hist(rvals.value, nbins)
bins, hist =hist_vals[1], hist_vals[0]

dist = np.exp(-bins/rscale.value)
dist = dist * np.max(hist) / np.max(dist)
plt.plot(bins, dist)
plt.xlim(0.0,600.0)
plt.savefig('dist.png')
plt.close()
