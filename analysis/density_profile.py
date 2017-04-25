import yt
import cgs as cgs
import dwarf as dw
import numpy as np
import matplotlib.pyplot as plt

# simulation directory
work_dir = '/home/emerick/Research/dwarfs/flash_runs/leo_T/LT_n075_v2_nh4_radtunnel/'
base_name = 'dwarf_fullp_hdf5'

nbins      = 30
field_list         = ['density','temperature','cell_mass']#['density']#,'temperature','pressure']
#field_ylims[field] = {'density' : , 'temperature': (1.0E3,1.0E6), 'pressure':} 

sim = dw.simulation(base_name, ds_dir = work_dir)              

xlim = (0.0, 600.0)
    
field_label = {'density' : r'Density (g cm$^{-3}$)', 'temperature': r'Temperature (K)', 'cell_mass' : r'Mass (M$_{\odot}$)'}
               
field_units = {'density' : 'g*cm**(-3)', 'temperature': 'K', 'cell_mass': 'Msun'}
xlabel = r'Radius (pc)'
for field in field_list:
    
    if field == 'cell_mass':
        accumulation=True ; weight_field = None
    else:
        accumulation=False; weight_field = 'cell_volume'
        
    
    dw.profile_1D(sim, 'radius', field, xlim*yt.units.pc, nbins, accumulation=accumulation, weight_field=weight_field) # compute the profile

    
    profile = sim.profiles['radius'][field]
    r       = sim.profile_bins['radius'].convert_to_units('pc').value
    t       = sim.profile_times#.convert_to_units('Myr').value
 
    # find max and mind values
    #minval = np.min(np.array(profile)) ; maxval = np.max(np.array(profile))
    
    maxval = np.min(profile[0].convert_to_units(field_units[field]).value)
    minval = np.max(profile[0].convert_to_units(field_units[field]).value)

    for p in profile:
        tempmin = np.min(p.convert_to_units(field_units[field]).value)
        tempmax = np.max(p.convert_to_units(field_units[field]).value)

        if minval > tempmin:
            minval = tempmin
        if maxval < tempmax:
            maxval = tempmax
            
    if field == 'temperature': 
        minval = 100.0 ; maxval = 6.0E7

    for i in np.arange(np.size(t)):
        
        plt.plot(r, profile[i].convert_to_units(field_units[field]).value, lw=2.5)
        plt.semilogy()
        
        plt.ylabel(field_label[field])
        plt.xlabel(xlabel)
        plt.ylim(minval, maxval)
        plt.xlim(xlim)
        plt.savefig(work_dir + 'profiles/' + field + '/' + field + '_radial_profile_%.1f.png'%(t[i].convert_to_units('Myr').value))
        plt.close()
    
        fname = work_dir + 'profiles/' + field + '/' + field + '_radial_profile_%.1f.dat'%(t[i].convert_to_units('Myr').value)
        np.savetxt(fname, np.transpose([r,profile[i].value]), fmt='%8.8E')
        
    
