import numpy as np
import yt

import os

import matplotlib.pyplot as plt
from   matplotlib    import rc

fsize = 17
rc('text', usetex=True)
rc('font', size=fsize)#, ftype=42)
line_width = 2.9


cfloor_only = True


work_dir = '/media/emerick/stardrive/emerick/simulations/FLASH/dwarfs/carina_final'
wdir_2   = '/home/emerick/Research/dwarfs/flash_runs/carina_final'
figure_name = 'carina_heating_cfloor.png'

files_dict = {'sne_9pc': wdir_2 + '/sne_9pc/carina_sne_9pc_mass.dat',
              'nosne_39pc_cfloor' : wdir_2 + '/nosne_38pc_cfloor/grav_bound_mass.dat',
              'sne_cfloor': wdir_2 + '/sne_9pc_cfloor/carina_sne_9pc_cfloor_mass.dat',
              'nosne_9pc' : wdir_2 + '/nosne_9pc/nosne_9pc_mass.dat',
              'nosne_cfloor' : wdir_2 + '/nosne_9pc_cfloor/carina_nosne_9pc_cfloor_mass.dat',
              'sne_19pc_cfloor': wdir_2 +'/sne_19pc_cfloor/grav_bound_mass.dat',
              'sne_78pc_cfloor': wdir_2 + '/sne_78pc_cfloor/grav_bound_mass.dat',
              'sne_39pc_cfloor': wdir_2 + '/sne_39pc_cfloor/grav_bound_mass.dat'}
              
              
fig, ax = plt.subplots(1,1) 
colors = {'sne_9pc' : 'black', 'nosne_9pc' : 'blue',
          'sne_cfloor' : 'black', 'nosne_cfloor' : 'blue', 
          'sne_19pc_cfloor':'green', 'sne_39pc_cfloor': 'red', 'sne_78pc_cfloor': 'orange',
          'nosne_39pc_cfloor' :'red'}

lstyle = {'sne_9pc' : '-', 'nosne_9pc' : '-',
          'sne_cfloor' : '--', 'nosne_cfloor' : '--', 'sne_19pc_cfloor':'--',
          'sne_39pc_cfloor' : '--', 'sne_78pc_cfloor': '--','nosne_39pc_cfloor':':'}

labels = {'sne_9pc' : 'H/C + SNe', 'nosne_9pc' : 'H/C',
          'sne_cfloor' : 'No H + SNe 9pc',
          'nosne_cfloor' : 'No H, No SNe', 'sne_19pc_cfloor' : 'No H + SNe 19pc',
          'sne_39pc_cfloor' : 'No H + SNe 39pc', 'sne_78pc_cfloor': 'No H + SNe 78pc',
          'nosne_39pc_cfloor' : 'No H - No SNe 39pc'}

        
xlabel = r'Time (Myr)'
ylabel = r'M$_{\rm{cold,gas}}$ / M$_{\rm{o}}$'

xlim = (0.0, 900.0)
ylim = (0.0,    1.1)             
    
if cfloor_only:
    files_dict.pop('sne_9pc', None)
    files_dict.pop('nosne_9pc', None)
    files_dict.pop('nosne_cfloor', None)
    
for name in files_dict:

    file_name = files_dict[name]
    data = np.genfromtxt(file_name, names=True)  
    
    M = data['m'] / data['m'][0]
    t = data['t']
    
    ax.plot(t, M, lw = line_width, color = colors[name], ls = lstyle[name], label = labels[name])
    
    
data = np.genfromtxt('./CarMedMidMass.dat', names = True)
ax.plot(data['t']*10.0, data['m'] / data['m'][0], color = 'red', lw = line_width,
                        ls='-', label = 'G13')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.legend(loc='upper right')                        
plt.tight_layout()
fig.savefig(figure_name)                        

             
