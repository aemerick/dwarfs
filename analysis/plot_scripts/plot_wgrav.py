import matplotlib.pyplot as plt
import numpy as np

from   matplotlib    import rc

import dwarf as dw

fsize = 17
rc('text', usetex=True)
rc('font', size=fsize)#, ftype=42)
line_width = 2.9

home = '/home/emerick/Research/dwarfs/flash_runs/leo_T/'

fig, ax = plt.subplots(1,2)

NFW_v2    = np.genfromtxt(home + 'SN_LT_n150_v2_nh4/0000_cfloor_global/grav_bound_mass.dat',names=True)
wgrav8_v2 = np.genfromtxt(home + 'SN_LT_n150_v2_nh4/0000_cfloor_wgrav/grav_bound_mass.dat',names=True)
wgrav9_v2 = np.genfromtxt(home + 'SN_LT_n150_v2_nh4/0000_cfloor_wgrav_2/grav_bound_mass.dat',names=True)

NFW_v4    = np.genfromtxt(home + 'SN_LT_n150_v4_nh4/0000_cfloor_global/grav_bound_mass.dat',names=True)
wgrav8_v4 = np.genfromtxt(home + 'SN_LT_n150_v4_nh4/0000_cfloor_wgrav/grav_bound_mass.dat', names=True)
wgrav9_v4 = np.genfromtxt(home + 'SN_LT_n150_v4_nh4/0000_cfloor_wgrav_2/grav_bound_mass.dat', names=True)

wgrav8_v2_tidal = np.genfromtxt(home + 'SN_LT_n150_v2_nh4/cfloor_tidal/grav_bound_mass.dat', names=True)


# 0 -> v = 2 plots
ax[0].plot(   NFW_v2['t'],       NFW_v2['m']/NFW_v2['m'][0], label = r'No Tidal Stripping', color = 'black', lw = 3, ls ='-')
ax[0].plot(wgrav9_v2['t'], wgrav9_v2['m']/wgrav9_v2['m'][0], label = r'0.9$\times \Phi$ - r$_{\rm{p}} \sim$ 100 kpc',
                                                             color = 'black', lw = 3, ls = '--')
ax[0].plot(wgrav8_v2['t'], wgrav8_v2['m']/wgrav8_v2['m'][0], label = r'0.8$\times \Phi$ - r$_{\rm{p}} \sim$ 30 kpc',
                                                             color = 'black', lw = 3, ls = ':')

#ax[0].plot(wgrav8_v2_tidal['t'], wgrav8_v2_tidal['m']/wgrav8_v2_tidal['m'][0], label = r'r$_{p}$ = 30 kpc - Tidal', color = 'green', lw = 3, ls ='-')

# ax 1 -> v = 400 plots
ax[1].plot(   NFW_v4['t'],       NFW_v4['m']/NFW_v4['m'][0], label = r'No Tidal Stripping', color = 'black', lw = 3, ls ='-')
ax[1].plot(wgrav9_v4['t'], wgrav9_v4['m']/wgrav9_v4['m'][0], label = r'0.9$\times \Phi$ - r$_{\rm{p}} \sim$ 100 kpc',
                                                             color = 'black', lw = 3, ls = '--')
ax[1].plot(wgrav8_v4['t'], wgrav8_v4['m']/wgrav8_v4['m'][0], label = r'0.8$\times \Phi$ - r$_{\rm{p}} \sim$ 30 kpc',
                                                             color = 'black', lw = 3, ls = ':')

##
print 'v2'
print 'NFW', dw.predict_stripping_time(NFW_v2['t'], NFW_v2['m'])
print '30 kpc', dw.predict_stripping_time(wgrav8_v2['t'], wgrav8_v2['m'])
print '100 kpc', dw.predict_stripping_time(wgrav9_v2['t'], wgrav9_v2['m'])

print 'v4'
print 'NFW', dw.predict_stripping_time(NFW_v4['t'], NFW_v4['m'])
print '30 kpc', dw.predict_stripping_time(wgrav8_v4['t'], wgrav8_v4['m']) 
print '100 kpc', dw.predict_stripping_time(wgrav9_v4['t'], wgrav9_v4['m'])

# set limits and labels
for a in ax:
    a.set_xlim(0.0, 2000.0)
    a.set_ylim(0.0, 1.1)
    a.set_ylabel(r'M(t) / M$_{\rm{o}}$')
    a.set_xlabel(r'Time (Myr)')
    a.minorticks_on()

    
ax[0].legend(loc='lower left')

xtext = 1300
ytext = 1.0
ax[0].annotate(r'v = 200 km s$^{-1}$', xy=(xtext, ytext), xytext=(xtext,ytext))
ax[1].annotate(r'v = 400 km s$^{-1}$', xy=(xtext, ytext), xytext=(xtext,ytext))

fig.set_size_inches(12,6)
plt.tight_layout()
fig.savefig('tidal_stripping.png')
