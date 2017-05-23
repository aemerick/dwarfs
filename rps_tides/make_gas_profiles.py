import profiles
import numpy as np
import yt
import deepdish as dd

def make_gas_profiles(ds = [], outname = 'gas_profile_evolution.h5'):

    all_data_dict = {}

    all_data_dict['t'] = np.zeros(len(ds))
    all_data_dict['surface_density'] = [None]*len(ds)
    all_data_dict['column_density']  = [None]*len(ds)

    for d in ds:
        com = profiles.center_of_mass(d)
        r, sigma, N = profiles.generate_gas_profile(d, d.all_data(), com = com)


        all_data_dict['r'] = r
        all_data_dict['t'][i] = d.current_time.convert_to_units('Myr').value
        all_data_dict['surface_density'][i] = sigma
        all_data_dict['column_density'][i]  = N


    dd.io.save(outname, all_data_dict)

    return all_data_dict

if __name__ == '__main__':

    dsnames = glob.glob('DD????/DD????')
    dsnames = np.sort(dsnames)

    make_gas_profiles(dsnames)
        
