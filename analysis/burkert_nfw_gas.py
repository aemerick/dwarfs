import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

fsize = 17
rc('text', usetex=True)
rc('font', size=fsize)#, ftype=42)
line_width = 2.9

# put rc here
#
#

mass_filename = 'grav_bound_mass.dat'
output_name = 'NFW_burkert_mass.png'


#mass_filename = 'coontained_mass.dat'
#output_name   = 'NFW_burkert_contained_mass.png'

# names
wdir = '/home/aemerick/Research/dwarfs/flash_runs/leo_T/'
filepaths = {'NFW_v2' : 'SN_LT_n150_v2_nh4/0000_cfloor_global/',
             'NFW_v4' : 'SN_LT_n150_v4_nh4/0000_cfloor_global/',
             'Burkert_v2' : 'Burkert_SN_LT_Mn150_v2_nh4/0000_cfloor_global/',
             'Burkert_v4' : 'Burkert_SN_LT_Mn150_v4_nh4/0000_cfloor_global/'}

for f in filepaths:
    filepaths[f] = wdir + filepaths[f] + mass_filename




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

xpos = 1300.0
ypos = 1.0

ax[0].annotate(r'v = 200 km s$^{-1}$', xy=(xpos,ypos),xytext=(xpos,ypos),
                                                    textcoords='data')
ax[1].annotate(r'v = 400 km s$^{-1}$', xy=(xpos,ypos),xytext=(xpos,ypos),
                                                    textcoords='data')
for a in ax:
    a.set_xlabel(r'Time (Myr)')
    a.set_ylabel(r'M(t) / M$_{\rm{o}}$')

    a.set_ylim(0.0,1.1)
    a.minorticks_on()
    a.set_xlim(0.0,2000.0)
    a.legend(loc='lower left')


    
fig.set_size_inches(12,6)
plt.tight_layout()
fig.savefig(output_name)
