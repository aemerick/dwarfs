import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
fsize = 15
rc('text')
rc('font', size=fsize)#, ftype=42)
line_width = 2.5

# put rc here
#
#


# names
wdir = '/home/emerick/Research/dwarfs/flash_runs/leo_T/'
filepaths = {'NFW_v2' : 'SN_LT_n150_v2_nh4/0000_cfloor_global/',
             'NFW_v4' : 'SN_LT_n150_v4_nh4/0000_cfloor_global/',
             'Burkert_v2' : 'Burkert_SN_LT_Mn150_v2_nh4/0000_cfloor_global/',
             'Burkert_v4' : 'Burkert_SN_LT_Mn150_v4_nh4/0000_cfloor_global/'}

for f in filepaths:
    filepaths[f] = wdir + filepaths[f] + 'grav_bound_mass.dat'




fig, ax = plt.subplots(1,2)

for f in ['NFW_v2','NFW_v4','Burkert_v2','Burkert_v4']:

    data = np.genfromtxt(filepaths[f], names=True)
    m    = data['m']
    t    = data['t']

    if 'v2' in f:
        ax_num = 0
    else:
        ax_num = 1

    if 'NFW' in f:
        color = 'black'; ls = '-'; label = 'NFW / Cusp'
    else:
        color = 'black'; ls = '--'; label = 'Burkert / Core'

    ax[ax_num].plot(t, m/m[0], lw = line_width, color = color, ls = ls, label = label)

ax[0].annotate(r'v$_{\rm{gal}}$ = 200 km s$^{-1}$', xy=(1000,1),xytext=(1000,1),
                                                    textcoords='data')
ax[1].annotate(r'v$_{\rm{gal}}$ = 400 km s$^{-1}$', xy=(1000,1),xytext=(1000,1),
                                                    textcoords='data')
for a in ax:
    a.set_xlabel(r'Time (Myr)')
    a.set_ylabel(r'M / M$_{o}$')

    a.set_ylim(0.0,1.1)

    a.set_xlim(0.0,2000.0)
    a.legend(loc='lower left')


    
fig.set_size_inches(12,6)
fig.savefig('NFW_burkert_mass.png')
