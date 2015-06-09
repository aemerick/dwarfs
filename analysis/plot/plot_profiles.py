import numpy as np
import matplotlib.pyplot as plt
import yt
import numpy as np

# some of my stuff:
import dwarf as dw # dwarf analysis code
import cgs as cgs  # cgs constants and conversions
import copy #?


ds_name = 'dwarf_fullp_hdf5'
ds_dir  = './../../flash_runs/carina_final/nosne_nowind_8pc/'
ds_selection = [0,50,100,150,200,221]

x_field = 'radius'    ; x_units = 'pc'
y_field = 'density'   ; y_units = 'g/cm**3'
time_units = 'Myr'

normalize = True
outname = 'nosne_nowind_' + x_field + "_" + y_field + '_profiles.png'

nbin = 20

# code to plot profiles of some simulation run:
sim = dw.simulation(ds_name, ds_dir=ds_dir,reload_times=False)

# make the profile
dw.profile_1D(sim, x_field, y_field, nbin, ds_selection = ds_selection)

y_list = copy.deepcopy(sim.profiles[x_field][y_field])
x_list = copy.deepcopy(sim.profile_bins[x_field])
t_list = copy.deepcopy(sim.profile_times)

#print x_list
#print y_list

#for i in np.arange(len(ds_selection)):
y_list = [x for x in y_list if x is not None]
t_list = [x for x in t_list if x is not None]


for i in np.arange(len(ds_selection)):
    x_list[i] = x_list[i].convert_to_units(x_units).value
    y_list[i] = y_list[i].convert_to_units(y_units).value
    t_list[i] = t_list[i].convert_to_units(time_units).value


#fig = plt.figure(figsize=[6,6])
#ax1 = fig.add_sublot(111)#

if normalize:
    normalization = y_list[0]
else:
    normalization = 1.0

for i in np.arange(len(ds_selection)):
    plt.plot(x_list[i], y_list[i] / normalization, lw=2, label="%.1f Myr"%(t_list[i]))


plt.xlabel(x_field + '('+x_units+')')
plt.ylabel(y_field + '(' + y_units +')')
plt.legend(loc='best')

plt.savefig(outname)

