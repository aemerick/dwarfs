import numpy as np
import matplotlib.pyplot as plt
import itertools

def flip(items, ncol):
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])

fpath = '/home/emerick/Research/dwarfs/flash_runs/'
lw = 2.0
ax = plt.subplot(111)

filenames = [fpath + "carina_adiabatic/carina_adiabatic_mass_dt5Myr.dat",
             fpath + "carina_final/adiabatic_31pc/carina_adiabatic_31pc_0-1000_dt5.dat",
             fpath + "carina_final/adiabatic_gatto_res/carina_adiabatic_gatto_res_0-900_dt5.dat",
             # no supernova:
             fpath + "carina_final/lowres_nosne/carina_nosne_bound_mass_0-750Myr_dt10.dat",
             fpath + "carina_final/nosne_gatto_res/carina_nosne_gatto_res_0-1000_dt5.dat",
             fpath + "carina_final/nosne_31pc/carina_nosne_31pc_0-1000_dt5.dat",
             # supernova :
             fpath + "carina_final/lowres/carina_nosne_mass_dt5Myr.dat",
             fpath + "carina_final/sne_31pc/carina_sne_31pc_0-512_dt2.dat",
             fpath + "carina_final/sne_gatto_res/carina_sne_gatto_res_0-1000_dt5.dat",
             # gatto
             "/home/emerick/Research/dwarfs/analysis/plot/CarMedMidMass.dat"]
             
             
labels = ['Adiabatic - 4.9pc', 'Adiabatic - 19pc', 'Adiabatic - 39 pc', 
          'No SNe - 4.9 pc', 'No SNe - 19 pc', 'No SNe - 39pc',
          'SNe - 4.9 pc', 'SNe - 19pc', 'SNe - 39pc',
          'Gatto 2D - 37 pc']

          
colors = ['blue', 'blue', 'blue', 'green', 'green', 'green', 'black', 'black', 'black', 'red']    
lines  = ['-', '--', '-.', '-', '--', '-.', '-', '--', '-.', '-']

data_sets = np.array([None]*len(labels))
for i in np.arange(len(labels)):
    data_sets[i] = np.genfromtxt(filenames[i], names=True)

# adjust gatto's data set
data_sets[-1]['t'] = data_sets[-1]['t'] * 10.0
data_sets[-1]['m'] = data_sets[-1]['m'] * 1.0E4
    
for i in np.arange(len(labels)):
    ax.plot(data_sets[i]['t'], data_sets[i]['m']/data_sets[i]['m'][0],
             label=labels[i], color=colors[i], ls=lines[i], lw=lw)
             

ax.set_xlabel('Time (Myr)')
ax.set_ylabel(r'M/M$_{o}$')
ax.set_ylim(0.0,1.35)
ax.set_xlim(0.0,1000.0)

hand, lab = ax.get_legend_handles_labels()
plt.legend(flip(hand, 2), flip(lab, 2), loc='best', ncol=2)
plt.savefig('carina_resolution_test.png')
plt.close()

