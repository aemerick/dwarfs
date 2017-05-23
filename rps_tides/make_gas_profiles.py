import profiles
import numpy as np
import yt
import deepdish as dd
import glob

def make_gas_profiles(dslist = [], outname = 'gas_profile_evolution.h5'):

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

if __name__ == '__main__':

    dsnames = glob.glob('DD????/DD????')
    dsnames = np.sort(dsnames)

    make_gas_profiles(dsnames)
