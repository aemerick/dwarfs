import numpy as np
import yt
import dwarf as dw
from initial_conditions import ic_list as icl
import os

import matplotlib.pyplot as plt
from   matplotlib    import rc

fsize = 15
rc('text', usetex=True)
rc('font', size=fsize)#, ftype=42)
line_width = 2.5


work_dir = '/home/emerick/Research/dwarfs/flash_runs/leo_T'
base_name = 'dwarf_fullp_hdf5'
figure_name = 'LT_adiabatic_mass_grav.png'

mass_file_ext = 'grav_bound_mass.dat'

# load known subdirectories
sub_dirs = next(os.walk(work_dir))[1]

sim_names = [s for s in sub_dirs if 'LT' in s and ((not 'radbox' in s) and (not 'SN' in s) and (not 'sn' in s) and (not 'radtunnel' in s))]
sim_dirs  = [work_dir + '/' + s for s in sim_names]

mass_files = [s + '/' + mass_file_ext for s,name in zip(sim_dirs,sim_names)]

# set up the panel plot
#
# make it a 3 row and 2 column plot
#
#
#
#
fig, ax = plt.subplots(3,2) 

n020_loc = 0 ; v2_loc = 0 ; nh3_color = 'blue'
n075_loc = 1 ; v4_loc = 1 ; nh4_color = 'black'
n150_loc = 2 ;              nh5_color = 'red'

nh3_label = r'n$_{\rm{halo}}$ = 10$^{%1i}$ cm$^{-3}$'%(-3)
nh4_label = r'n$_{\rm{halo}}$ = 10$^{%1i}$ cm$^{-3}$'%(-4)
nh5_label = r'n$_{\rm{halo}}$ = 10$^{%1i}$ cm$^{-3}$'%(-5)


xlabel = r'Time (Myr)'
ylabel = r'M$_{\rm{cold,gas}}$ / M$_{\rm{o}}$'

xlim = (0.0, 2000.0)
ylim = (0.0,    1.0)

for i in np.arange(len(sim_names)):
    name = sim_names[i]
    sim_dir = sim_dirs[i]

    if 'n020' in name:
        plot_row = n020_loc ; model_label = r'n$_{\rm{o}}$ = %.2f'%(0.02)
    elif 'n075' in name:
        plot_row = n075_loc ; model_label = r'n$_{\rm{o}}$ = %.2f'%(0.75)
    else:
        plot_row = n150_loc;  model_label = r'n$_{\rm{o}}$ = %.2f'%(1.50)
    if 'v2' in name:
        plot_column = v2_loc
    else:
        plot_column = v4_loc
    if 'nh3' in name:
        color = nh3_color
        label = nh3_label 
    elif 'nh4' in name:
        color = nh4_color
        label = nh4_label #r'n$_{\rm{halo}}$ = 10$^{%1i}$ cm$^{-3}$'%(-4)

    else:
        color = nh5_color
        label = nh5_label #r'n$_{\rm{halo}}$ = 10$^{%1i}$ cm$^{-3}$'%(-5)
    

    plot_loc = (plot_row, plot_column)

    mass_file = mass_files[i]

    if os.path.isfile(mass_file):

        try:
            data = np.genfromtxt(mass_file,names=True)
            ax[plot_loc].plot(data['t'] , data['m']/data['m'][0], color = color,
                                          lw = line_width, ls = '-', label=label)

            axis = ax[plot_loc]
            axis.annotate(model_label,xy=(1400,0.8),xytext=(1400,0.8),textcoords='data')

          #  print name, data['t'][ np.where(data['m'][data['t']<100] == np.min(data['m'][data['t'] < 100]))]
            
            
            t_strip = dw.predict_stripping_time(data['t'], data['m'], t_fit = 100.0)

            print name + ' stipping time %.3f Gyr'%(t_strip/1000.0)
            
            
            
        except:
            print 'nothing'        
    else:
        print 'Mass evolution file does not exist for ' + name

    ax[plot_loc].set_xlim(xlim)
    ax[plot_loc].set_ylim(ylim)
    ax[plot_loc].set_xlabel(xlabel)
    ax[plot_loc].set_ylabel(ylabel)

ax[(0,v2_loc)].set_title(r'v$_{\rm{gal}}$ = %3i km s$^{-1}$'%(200.0))
ax[(0,v4_loc)].set_title(r'v$_{\rm{gal}}$ = %3i km s$^{-1}$'%(400.0))



#ax[0,0].plot(-1,-1,color=nh3_color,label=nh3_label, ls='-', lw=line_width)
#ax[0,0].plot(-1,-1,color=nh4_color,label=nh4_label, ls='-', lw=line_width)
#ax[0,0].plot(-1,-1,color=nh5_color,label=nh5_label, ls='-', lw=line_width)
#ax[(0,0)].legend(loc='best',fancybox=True)
#ax[0,0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
#          ncol=3, fancybox=True, shadow=True)
#fig.save()
plt.tight_layout()
fig.savefig(figure_name)
