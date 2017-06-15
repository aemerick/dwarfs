import matplotlib.pyplot as plt
import numpy as np
import cgs as cgs
#import dwarf_model as dw_model
from initial_conditions import ic_list as icl ;
from matplotlib import rc

from astropy import units as u

fsize = 15
rc('text', usetex=True)
rc('font', size=fsize)#, ftype=42)
line_width = 2.5
point_size = 30


# load all of the simulations from the initial conditions
LT_names = []
LT_dict  = {}
for ic in icl.ic_object_dict.keys():
    if ic.startswith('LT_n'):
        LT_names.append(ic)
        LT_dict[ic] = icl.ic_object_dict[ic]
#

# plot the points
nh3_color = 'blue' ; nh4_color = 'black'; nh5_color = 'red'

v2_ps = 'D' ; v4_ps = 's'


for name in LT_names:

    if 'nh3' in name:
        color = nh3_color
    elif 'nh4' in name:
        color = nh4_color
    else:
        color = nh5_color

    if 'v2' in name:
        marker = v2_ps
    else:
        marker = v4_ps    

    P_RPS = LT_dict[name].ic['n_halo'] * cgs.mp * LT_dict[name].ic['mu_halo'] * LT_dict[name].ic['v_halo']**2.0

    
    
    plt.scatter( LT_dict[name].ic['n_o'] , P_RPS, marker = marker, s = point_size ,color= color)


# plot the alphas
alpha = [5.0, 3.0, np.pi/2.0, 1.0]
linestyle = ['-', '--', '-.', ':']
color = ['black','gray','purple']
#n_sample = np.linspace(0.0,2.0,50)
#for i in np.arange(len(alpha)):
# 
#    j = 0
#    for name in ['LT_n075_v2_nh5']:
#        M_o = LT_dict[name].ic['M_HI'] + LT_dict[name].ic['M_DM']
#        r_o = LT_dict[name].ic['r_HI']
#
#        P_anal   = alpha[i] * cgs.G * n_sample * M_o / (r_o) * 1.31 * cgs.mp
#        plt.plot(n_sample, P_anal, ls = linestyle[i], lw = 2.5, color = 'black')
#
#        j = j + 1

p_stable = np.genfromtxt('./../analytic_model/p_ram_stable.dat',names=True)



plt.plot(p_stable['n'], p_stable['10'], ls = '-', lw = 2.5, color = 'black')
plt.plot(p_stable['n'], p_stable['25'], ls = '--', lw = 2.5, color = 'black')
plt.plot(p_stable['n'], p_stable['50'], ls = '-.', lw = 2.5, color = 'black')
plt.plot(p_stable['n'], p_stable['100'], ls = ':', lw = 2.5,  color = 'black')

#plot Slater Bell 2014 lines
SB_50 = 10.0**(-12.8)
SB_90 = 10.0**(-14.8)

plt.plot(plt.xlim(), 2*[SB_50], ls='-', lw=2.5, color='green')
plt.plot(plt.xlim(), 2*[SB_90], ls='--', lw=2.5, color='green')

plt.ylim(1.0E-18,5.0E-12)
plt.xlim(0.0,2.0)
plt.semilogy()
plt.ylabel(r'RPS Pressure (dyne cm$^{-2}$)')
plt.xlabel(r'Dwarf Galaxy - Central Gas Density (cm$^{-3}$)')


# do legend things
plt.scatter(-1,-1, color = nh3_color, marker='o', s=point_size, 
                  label = r'n$_{\rm{halo}}$ = 10$^{-3}$ cm$^{-3}$')
plt.scatter(-1,-1, color = nh4_color, marker='o', s=point_size, 
                  label = r'n$_{\rm{halo}}$ = 10$^{-4}$ cm$^{-3}$')
plt.scatter(-1,-1, color = nh5_color, marker='o', s=point_size, 
                  label = r'n$_{\rm{halo}}$ = 10$^{-5}$ cm$^{-3}$')

plt.scatter(-1,-1, color ='black', marker = v2_ps, s= point_size,
                  label = r'v$_{\rm{galaxy}}$ = 200 km s$^{-1}$')
plt.scatter(-1,-1, color ='black', marker = v4_ps, s= point_size,
                  label = r'v$_{\rm{galaxy}}$ = 400 km s$^{-1}$')

plt.plot(-1,-1,color='black',ls='-',lw=2.5,label = r'Minimum P$_{\rm{ram}}$')

plt.legend(loc='best',fancybox=True,ncol=2)

plt.savefig('LT_RPS_condition.png')

