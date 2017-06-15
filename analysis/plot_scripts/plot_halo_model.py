import numpy as plt
from halo import galaxy as gal
import numpy as np
import cgs

from matplotlib import rc
import matplotlib.pyplot as plt
from scipy.optimize import brentq


line_width = 2.5 ; point_size = 35
fsize = 15
rc('text', usetex=True)
rc('font', size=fsize)#, ftype=42)
rc('legend', numpoints=1)
line_width = 2.5

def MB15(r,beta=0.5):
    #n_o * rc^(3beta)
    #r^(3beta)
    numerator = 1.35E-2
    denominator = r**(3.0*beta)
    return numerator / denominator

# sample r points
r = np.linspace(0.0, 250.0, 1000.0)* cgs.kpc
n_no_ambient         = gal.halo_gas_ndensity(r, ambient = 0.0)
n_ambient            = gal.halo_gas_ndensity(r)

# 
plt.plot(r/cgs.kpc, MB15(r/cgs.kpc), label='MB15',color='black',ls='-.',lw=line_width)


# plot MB13 lines
#plt.plot(r/cgs.kpc, n_no_ambient, label=r'MB13',  color = 'black', ls='-', lw=line_width)
plt.plot(r/cgs.kpc, n_ambient, label=r'MB13', color = 'black', ls='-',lw=line_width)
#plt.plot(r/cgs.kpc, n_ambient - n_no_ambient, label=r'MB13 Ambient', color='black',ls='--',lw=line_width)
# Grcevich & putman 2009
#plt.plot([1.0,70.0], [2.0E-4,2.0E-4], label='GP09',
#                                color='green', lw=line_width, ls='-')
#plt.plot([1.0,70.0],[3.0E-4,3.0E-4],color='green', lw=line_width, ls='-')

# lines from Grcevich and Putman 2009 --- just points, no error bars yet
n_GP09 = [8.5E-5, 2.1E-4, 2.7E-4, 3.1E-4]
r_GP09 = [20    ,     40,     68,    118]
plt.scatter(r_GP09, n_GP09, label = 'GP09', color= 'blue', s = point_size)

low_r_err = [17.0, 30.0, 35.0, 62.0]
up_r_err  = [43.0, 36.0, 15.0, 26.0]
low_n_err = [0.3E-4,  0.08E-4, 2.26E-4, 2.12E-4]
up_n_err  = [3.05E-4, 5.1E-4, 1.2E-4, 1.5E-4]

plt.errorbar(r_GP09, n_GP09, xerr = [low_r_err, up_r_err], yerr = [low_n_err, up_n_err], lw=line_width,color='blue')
    
# plot Gatto points
# 
r_gatto = [73.5, 64.7] 
lower_r_error = [13.7, 13.5]
upper_r_error = [16.7,17.01]             
n_gatto = [3.15E-4, 2.55E-4]
n_err_gatto = [1.85E-4, 1.05E-4]
plt.errorbar(r_gatto, n_gatto, xerr=[lower_r_error, upper_r_error], yerr=n_err_gatto,color='orange', lw=line_width)
plt.scatter(r_gatto,n_gatto,label='G13',color='orange',s=point_size, marker='s')

#
# plot Salem 2015
#
plt.errorbar(48.2, 1.1E-4, xerr=0.5, yerr=0.45E-4, color='black' ,lw=line_width)
plt.scatter(48.2, 1.1E-4, label='S15', color='black',marker='D', s=point_size)

# kaufmannn
#plt.plot(gal.kaufmann('01','r'), gal.kaufmann('01','density'), label = 'Kaufmann',
#                                           color = 'purple', lw=line_width, ls='-')

plt.plot(gal.kaufmann('02','r'), gal.kaufmann('02','density'), label = 'K09',
                                           color = 'purple', lw=line_width, ls='-')

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
#for r_value in r_sample:
 #   plt.plot( [r_value/cgs.kpc]*2 , plt.ylim(), color = 'black', ls=':', lw=line_width)


#plt.arrow(40, 2.1E-4, 40, 6.0E-3, fc='green', ec='green')
plt.minorticks_on()
plt.legend(loc='best',fancybox=True,ncol=2)
plt.xlim(np.min(r)/cgs.kpc,np.max(r)/cgs.kpc)
plt.ylim(2.0E-6,5.0E-3)
plt.xlabel(r'Radius (kpc)')
plt.ylabel(r'n$_{\rm{e}}$ (cm$^{-3}$)')
plt.savefig('hot_halo_models.png')
