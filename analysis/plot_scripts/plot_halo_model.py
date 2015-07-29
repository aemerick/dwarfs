import numpy as plt
from halo import galaxy as gal
import numpy as np
import cgs

from matplotlib import rc
import matplotlib.pyplot as plt
from scipy.optimize import brentq

line_width = 2.5
fsize = 15
rc('text', usetex=True)
rc('font', size=fsize)#, ftype=42)
line_width = 2.5


# sample r points
r = np.linspace(0.0, 250.0, 1000.0)* cgs.kpc
n_no_ambient         = gal.halo_gas_ndensity(r, ambient = 0.0)
n_ambient            = gal.halo_gas_ndensity(r)

# plot MB13 lines
plt.plot(r/cgs.kpc, n_no_ambient, label=r'MB13 - No Constant',  color = 'black', ls='--', lw=line_width)
plt.plot(r/cgs.kpc, n_ambient, label=r'MB13 - With Ambient', color = 'black', ls='-',lw=line_width)

# Grcevich & putman 2009
#plt.plot([1.0,70.0], [2.0E-4,2.0E-4], label='GP09',
#                                color='green', lw=line_width, ls='-')
#plt.plot([1.0,70.0],[3.0E-4,3.0E-4],color='green', lw=line_width, ls='-')

n_GP09 = [8.5E-5, 2.1E-4, 2.7E-4, 3.1E-4]
r_GP09 = [20    ,     40,     68,    118]
plt.plot(r_GP09, n_GP09, label = 'GP09', color= 'green', lw = line_width, ls='-')
#plt.arrow(40, 2.1E-4, 40, 6.0E-4, fc='green', ec='green')

# plot Gatto points

# kaufmannn
plt.plot(gal.kaufmann('01','r'), gal.kaufmann('01','density'), label = 'Kaufmann',
                                           color = 'purple', lw=line_width, ls='-')

plt.plot(gal.kaufmann('02','r'), gal.kaufmann('02','density'), label = 'Kaufmann Realistic',
                                           color = 'purple', lw=line_width, ls='--')

# plot radii corresponding to model

n_sample = [1.0E-3, 1.0E-4, 1.0E-5]
r_sample = [None]*3

i = 0
for n in n_sample:
    root_function = lambda x : n - gal.halo_gas_ndensity(x, ambient = 0.0)

    r_sample[i] = brentq(root_function, 0.0, 500.0*cgs.kpc)
    i = i  + 1

print n_sample, np.array(r_sample)/cgs.kpc
# plot them as verticle lines
plt.semilogy()
for r_value in r_sample:
    plt.plot( [r_value/cgs.kpc]*2 , plt.ylim(), color = 'blue', ls='-', lw=line_width)


#plt.arrow(40, 2.1E-4, 40, 6.0E-3, fc='green', ec='green')


plt.legend(loc='best',fancybox=True)

plt.xlabel(r'Radius (kpc)')
plt.ylabel(r'n$_{\rm{e}}$ (cm$^{-3}$)')
plt.savefig('hot_halo_models.png')
