import numpy as np
import yt

import os

import matplotlib.pyplot as plt
from   matplotlib    import rc

fsize = 17
rc('text', usetex=True)
rc('font', size=fsize)#, ftype=42)
line_width = 2.9

work_dir = '/media/stardrive/emerick/simulations/FLASH/dwarfs/carina_final'
wdir_2   = '/home/emerick/Research/dwarfs/flash_runs/carina_final'
figure_name = 'carina_resolution_study.png'

ad_files = {'38': work_dir + '/adiabatic_31pc/carina_adiabatic_31pc_0-1000_dt5.dat',
            '19': work_dir + '/adiabatic_gatto_res/carina_adiabatic_gatto_res_0-900_dt5.dat',
            '9' : wdir_2   + '/adiabatic_9pc/carina_adiabatic_9pc_mass.dat',
            '4' : wdir_2 + '/../carina_adiabatic/carina_adiabatic_mass_dt5Myr.dat'}

sne_files = {'38' : work_dir + '/sne_31pc/carina_sne_31pc_0-512_dt2.dat',
             '19' : work_dir + '/sne_gatto_res/carina_sne_gatto_res_0-1000_dt5.dat',
             '9'  : wdir_2   + '/sne_9pc/carina_sne_9pc_mass.dat',
             '4'  : work_dir + '/lowres/carina_nosne_mass_dt5Myr.dat'}
             
nosne_files = {'38' : work_dir + '/nosne_31pc/carina_nosne_31pc_0-1000_dt5.dat',
               '19' : work_dir + '/nosne_gatto_res/carina_nosne_gatto_res_0-1000_dt5.dat',
               '9'  : wdir_2   + '/nosne_9pc/nosne_9pc_mass.dat',
               '4'  : work_dir + '/lowres_nosne/carina_nosne_bound_mass_0-750Myr_dt10.dat'}
# set up the panel plot
#
# make it a 3 row and 2 column plot
#
#
#
#
fig, ax = plt.subplots(1,3) 

loc = {'ad':0,'nosne':1,'sne':2}
fdict = {'ad':ad_files,'nosne':nosne_files,'sne':sne_files}

colors = {'4' : 'black', '9': 'green', '19': 'blue', '38':'orange'}
lstyle = {'4' : '-'    , '9': '--',  '19':':', '38':'-.'}

labels = {'4' : '4.9 pc', '9' : '9.8 pc', '19' : '19.5 pc' , '38' : '39.1 pc'}

xlabel = r'Time (Myr)'
ylabel = r'M$_{\rm{cold,gas}}$ / M$_{\rm{o}}$'

xlim = (0.0, 900.0)
ylim = (0.0,    1.1)


for name in loc:
    pn = loc[name]

    for res in ['4','9','19','38']:
        file_name = fdict[name][res]
        data = np.genfromtxt(file_name, names=True)  
    
        M = data['m'] / data['m'][0]
        t = data['t']
    
        ax[pn].plot(t,M, lw = line_width, color = colors[res], ls = lstyle[res], label = labels[res])

    ax[pn].set_xlim(xlim)
    ax[pn].set_ylim(ylim)
    ax[pn].set_xlabel(xlabel)
    ax[pn].set_ylabel(ylabel)
    #ax[pn].legend(loc='upper right')
x = 35
y = 0.08

ax[loc['ad']].legend(loc='lower left') ; ax[loc['nosne']].legend(loc='upper right') ; ax[loc['sne']].legend(loc='upper right')

ax[loc['ad']].annotate(r'Adiabatic', xy=(x,y), xytext=(plt.xlim()[1]-300,plt.ylim()[1]-0.1), color='black')  
ax[loc['nosne']].annotate(r'Heating + Cooling', xy=(x,y), xytext=(x,y), color='black')     
ax[loc['sne']].annotate(r'H/C + SNe', xy=(x,y), xytext=(x,y), color='black')     

### plot cfloor lines 
#data = np.genfromtxt(wdir_2 + '/nosne_9pc_cfloor/carina_nosne_9pc_cfloor_mass.dat',names=True)
#ax[loc['nosne']].plot(data['t'], data['m']/data['m'][0], color = 'purple', lw=line_width)

### plot cfloor lines 
#data = np.genfromtxt(wdir_2 + '/sne_9pc_cfloor/carina_sne_9pc_cfloor_mass.dat',names=True)
#ax[loc['sne']].plot(data['t'], data['m']/data['m'][0], color = 'purple', lw=line_width)

fig.set_size_inches(14,6)
plt.tight_layout()
fig.savefig(figure_name)
