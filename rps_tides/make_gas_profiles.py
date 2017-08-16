import matplotlib.pyplot as plt
from matplotlib import rc
fsize = 17
rc('text', usetex=False)
rc('font', size=fsize)#, ftype=42)
line_width = 4
point_size = 30
viridis    = matplotlib.cm.get_cmap('viridis')

import profiles
import numpy as np
import yt
import deepdish as dd
import glob
import os


def make_dm_profiles(dslist = [], outname = 'dm_profile_evolution.h5', overwrite = False):

    if os.path.exists(outname):
        print "Output file exists - do not overwrite"
        return dd.io.load(outname)

    all_data_dict = {}

    all_data_dict['t'] = np.zeros(len(dslist))
    all_data_dict['cumulative_mass'] = [None]*len(dslist)
    all_data_dict['surface_density'] = [None]*len(dslist)

    i = 0
    for d in dslist:
        ds = yt.load(d)
        com = profiles.center_of_mass(ds)
        r, sigma, mass_r, M, rho = profiles.generate_dm_profile(ds, ds.all_data(), com=com)

        all_data_dict['r']    = r
        all_data_dict['r_sp'] = mass_r
        all_data_dict['t'][i] = ds.current_time.convert_to_units('Myr').value
        all_data_dict['surface_density'][i] = sigma
        all_data_dict['cumulative_mass'][i] = M
        all_data_dict['density'][i] = rho

        i = i + 1

    dd.io.save(outname, all_data_dict)

    return

def make_gas_profiles(dslist = [], outname = 'gas_profile_evolution.h5', overwrite = False):

    if os.path.exists(outname):
        print "Output file exists - do not overwrite"
        return dd.io.load(outname)

    all_data_dict = {}

    all_data_dict['t'] = np.zeros(len(dslist))
    all_data_dict['surface_density'] = [None]*len(dslist)
#    all_data_dict['column_density']  = [None]*len(dslist)

    i = 0
    for d in dslist:
        ds = yt.load(d)
        com = profiles.center_of_mass(ds)
        r, sigma = profiles.generate_gas_profile(ds, ds.all_data(), com = com)


        all_data_dict['r'] = r
        all_data_dict['t'][i] = ds.current_time.convert_to_units('Myr').value
        all_data_dict['surface_density'][i] = sigma

        i = i + 1
#        all_data_dict['column_density'][i]  = N


    dd.io.save(outname, all_data_dict)

    return all_data_dict

def plot_gas_profiles(data):

    fig, ax = plt.subplots()

    max_time = data['t'][-1]
    z        = data['t'] / max_time
    x        = (data['r'][1:] + data['r'][:-1])*0.5    

    for i in np.arange(np.size(data['t'])):
        ax.plot( x, data['surface_density'][i],
                     color = viridis(z[i]), lw = line_width)

    ax.set_ylabel(r'Gas Surface Density (M$_{\odot]$ pc$^{-2}$)')
    ax.set_xlabel(r'Radius (kpc)')

    fig.set_size_inches(8,8)
    plt.tight_layout()
    plt.minorticks_on()

    fig.savefig('gas_profile_evolution.png')

    return

def plot_dm_profiles(data):

    fig, ax = plt.subplots()

    max_time = data['t'][-1]
    z        = data['t'] / max_time
    x        = (data['r_sp'][1:] + data['r_sp'][:-1])*0.5

    for i in np.arange(np.size(data['t'])):
        ax.plot( x, data['density'][i],
                     color = viridis(z[i]), lw = line_width)

    ax.set_ylabel(r'Dark Matter Density (M$_{\odot]$ kpc$^{-3}$)')
    ax.set_xlabel(r'Radius (kpc)')

    fig.set_size_inches(8,8)
    plt.tight_layout()
    plt.minorticks_on()

    fig.savefig('dm_profile_evolution.png')

    return

if __name__ == '__main__':

    dsnames = glob.glob('DD????/DD????')
    dsnames = np.sort(dsnames)

    dm_profiles  = make_dm_profiles(dsnames)
    gas_profiles = make_gas_profiles(dsnames)

    plot_gas_profiles(gas_profiles)
    plot_dm_profiles(dm_profiles)
